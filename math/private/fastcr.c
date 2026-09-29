/* Fast front end for crlibm's exp_rd and log_rz (the functions MESA's crmath exp() and log() call).
 *
 * Both crlibm functions are correctly rounded (exp: toward -inf, log: toward zero), so the result for
 * every input is a unique double. Two layers, each returning exactly crlibm's bits:
 *  1. a per-thread 4-entry cache of recent (argument bits -> result). MESA repeats arguments a lot:
 *     auto_diff pow(x_ad, y) evaluates pow(x, y), pow(x, y-1), pow(x, y-2), pow(x, y-3), i.e. log(x)
 *     four times; Skye evaluates constant pows per call. RGB: ~39% of log and ~22% of exp calls hit.
 *  2. a double-double evaluation of exp/log with relative error far below 2^-63 (analysis <~2^-70,
 *     measured 2^-73 log / 2^-78 exp). The directed rounding is decided from the exact residual of
 *     the final double-double -> double step. If the residual is within 2^-63 of a rounding boundary
 *     (~0.1-1% of arguments, and every exact case), crlibm's own result is returned instead.
 *
 * Build: -O2 -ffp-contract=off (the error-free transformations must not be FMA-contracted; the
 * explicit FMAs come from target("fma") on the evaluation functions, used only if the CPU has FMA).
 * Production: linked into the star executable (math_init references savethesun_fastcr_link, which pulls
 * this object out of libmath.a) or LD_PRELOADed. The definitions of exp_rd / log_rz interpose the crlibm
 * ones called from libcrmath.so; crlibm's are reached with dlsym(RTLD_NEXT).
 * On by default; MESA_FAST_CRLIBM=0 makes every call forward to crlibm.
 * Test (-DFASTCR_TEST): defines fast_exp_rd / fast_log_rz, fallback calls crlibm directly, counts
 * fallbacks and exposes the last double-double approximation.
 */
#define _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>
#include "fastcr_tables.h"

#if defined(__x86_64__) || defined(__i386__)
#define TGT __attribute__((target("fma"), noinline))
#define HAVE_FMA_CPU() (__builtin_cpu_init(), __builtin_cpu_supports("fma"))
#else
#define TGT __attribute__((noinline))
#define HAVE_FMA_CPU() 1
#endif
#define FMA(a, b, c) __builtin_fma(a, b, c)   /* hardware FMA inside TGT functions; always exact */

void savethesun_fastcr_link(void) {}         /* referenced from math_init so the linker keeps this object */

static inline uint64_t asu(double x) { uint64_t u; memcpy(&u, &x, 8); return u; }
static inline double asd(uint64_t u) { double x; memcpy(&x, &u, 8); return x; }
static inline double fabs_(double x) { return asd(asu(x) & 0x7fffffffffffffffULL); }

#define EPS 0x1p-63

#ifdef FASTCR_TEST
double log_rz(double), exp_rd(double);
long fastcr_fb_log, fastcr_fb_exp;
double fastcr_dbg_h, fastcr_dbg_l;
#define FB_LOG(x) do { fastcr_fb_log++; return log_rz(x); } while (0)
#define FB_EXP(x) do { fastcr_fb_exp++; return exp_rd(x); } while (0)
#define REAL_LOG log_rz
#define REAL_EXP exp_rd
#define DBG(h, l) (fastcr_dbg_h = (h), fastcr_dbg_l = (l))
#define NAME_LOG fast_log_rz
#define NAME_EXP fast_exp_rd
#define ENABLED 1
#else
typedef double (*fn_t)(double);
static fn_t real_log_rz, real_exp_rd;
static int fastcr_on;
static double missing(double x) { (void)x; fprintf(stderr, "fastcr: crlibm log_rz/exp_rd not found\n"); abort(); }
__attribute__((constructor)) static void fastcr_init(void) {
  /* also runs in processes without crlibm (e.g. when preloaded): only complain if actually called */
  real_log_rz = (fn_t)dlsym(RTLD_NEXT, "log_rz");
  real_exp_rd = (fn_t)dlsym(RTLD_NEXT, "exp_rd");
  if (!real_log_rz || !real_exp_rd) { real_log_rz = real_exp_rd = missing; return; }
  const char *e = getenv("MESA_FAST_CRLIBM");
  if (e && e[0] == '0') return;
  /* silent when on by default: MESA's tests diff program output (math/test/ck keys on its first line) */
  if (!HAVE_FMA_CPU()) {
    if (e) { printf(" savethesun: MESA_FAST_CRLIBM ignored (CPU has no FMA)\n"); fflush(stdout); }
    return;
  }
  fastcr_on = 1;
  if (e) { printf(" savethesun: fast correctly-rounded exp/log front end (crlibm fallback)\n"); fflush(stdout); }
}
#define FB_LOG(x) return real_log_rz(x)
#define FB_EXP(x) return real_exp_rd(x)
#define REAL_LOG real_log_rz
#define REAL_EXP real_exp_rd
#define DBG(h, l) ((void)0)
#define NAME_LOG log_rz
#define NAME_EXP exp_rd
#define ENABLED fastcr_on
#endif

/* ---- per-thread cache of the last 4 distinct arguments (keys are argument bits) */
typedef struct { uint64_t k[4]; double v[4]; unsigned i; } cache4;
#define TLS __thread __attribute__((tls_model("initial-exec")))   /* no __tls_get_addr calls */
static TLS cache4 c_log = {{0x3ff0000000000000ULL, 0x3ff0000000000000ULL, 0x3ff0000000000000ULL, 0x3ff0000000000000ULL},
                                {0.0, 0.0, 0.0, 0.0}, 0};          /* log_rz(1) = +0 */
static TLS cache4 c_exp = {{0, 0, 0, 0}, {1.0, 1.0, 1.0, 1.0}, 0};  /* exp_rd(+0) = 1 */
static inline int cget(const cache4 *c, uint64_t b, double *r) {
  int q = c->k[0] == b ? 0 : c->k[1] == b ? 1 : c->k[2] == b ? 2 : c->k[3] == b ? 3 : -1;
  if (q < 0) return 0;
  *r = c->v[q]; return 1;
}
static inline double cput(cache4 *c, uint64_t b, double r) {
  c->i = (c->i + 1) & 3; c->k[c->i] = b; c->v[c->i] = r; return r;
}

/* ---- exp(x) rounded toward -inf */
static inline __attribute__((always_inline)) double exp_core(double x) {
  uint64_t ax = asu(x) & 0x7fffffffffffffffULL;
  if (ax == 0) return 1.0;                      /* exact case (crlibm: 1 for +-0) */
  if (ax > 0x4085e00000000000ULL) FB_EXP(x);    /* |x| > 700 (keeps the result normal), inf, nan */
  /* x = k ln2/128 + r, |r| <~ ln2/256; r = rh + rl exactly up to ~2^-78 */
  const double shift = 0x1.8p52;
  double kd = x * EXP_INVL + shift;
  int64_t ki = (int64_t)(asu(kd) - asu(shift));
  kd -= shift;
  double a = x - kd * EXP_LH;                   /* exact: k*LH has <= 53 bits, Sterbenz */
  double b = kd * EXP_LL;
  double rh = a - b;
  double bb = rh - a, rl = (a - (rh - bb)) + (-b - bb);   /* TwoSum(a, -b) */
  /* expm1(r) = rh + rh^2/2 + [rl + rh*rl + (rh^2 err)/2 + rh^3 P(rh)] */
  double sh = rh * rh, sl = FMA(rh, rh, -sh);
  /* P(r) = 1/6 + r/24 + r^2/120 + r^3/720 + r^4/5040, Estrin with FMA (short dependency chain) */
  double q01 = FMA(rh, 0x1.5555555555555p-5, 0x1.5555555555555p-3);
  double q23 = FMA(rh, 0x1.6c16c16c16c17p-10, 0x1.1111111111111p-7);
  double p3 = (sh * rh) * FMA(sh, FMA(sh, 0x1.a01a01a01a01ap-13, q23), q01);
  double hs = 0.5 * sh;
  double eh = rh + hs, el = hs - (eh - rh);     /* Fast2Sum, |rh| >= |hs| */
  double lo = ((el + rl) + (rh * rl + 0.5 * sl)) + p3;
  /* 2^(j/128) (1 + E) */
  int j = (int)(ki & 127); int64_t m = ki >> 7;
  double Th = exp_T[2 * j], Tl = exp_T[2 * j + 1];
  double ph = Th * eh, pl = FMA(Th, eh, -ph);
  double h = Th + ph, hl = ph - (h - Th);       /* Fast2Sum, Th >= 1 > |ph| */
  double low = ((hl + pl) + (Tl + Tl * eh)) + Th * lo;
  double u = h + low, t = low - (u - h);        /* h + low = u + t exactly */
  DBG(h, low);
  if (!(fabs_(t) > EPS * u)) FB_EXP(x);
  if (t < 0) u = asd(asu(u) - 1);               /* round down: predecessor (u > 0) */
  return u * asd((uint64_t)(1023 + m) << 52);   /* exact scaling, normal range */
}

/* ---- log(x) rounded toward zero */
static inline __attribute__((always_inline)) double log_core(double x) {
  uint64_t ix = asu(x);
  if (ix == 0x3ff0000000000000ULL) return 0.0;  /* exact case (crlibm: +0) */
  if (ix - 0x0010000000000000ULL >= 0x7ff0000000000000ULL - 0x0010000000000000ULL) FB_LOG(x); /* <=0, subnormal, inf, nan */
  /* x = 2^k z, z in [0.6875, 1.375), bucket i; r = z*invc - 1 = rh + rl exactly */
  uint64_t tmp = ix - 0x3fe6000000000000ULL;
  int i = (int)((tmp >> 45) & 127);
  int64_t k = (int64_t)tmp >> 52;
  double z = asd(ix - (tmp & (0xfffULL << 52)));
  double invc = log_T[3 * i], Lh = log_T[3 * i + 1], Ll = log_T[3 * i + 2];
  double ph = z * invc, pl = FMA(z, invc, -ph);
  double rh = ph - 1.0, rl = pl;                /* ph - 1 exact (Sterbenz) */
  double kd = (double)k;
  /* log1p(r) = r - r^2/2 + r^3/3 + D, D = sum_{n=4..10} (-1)^(n+1) r^n / n */
  double sh = rh * rh, sl = FMA(rh, rh, -sh);
  double t3h = sh * rh, t3l = FMA(sh, rh, -t3h) + sl * rh;
  double ch = t3h * THIRD_H, cl = (FMA(t3h, THIRD_H, -ch) + t3h * THIRD_L) + t3l * THIRD_H;
  double rr = rh + rl, r2 = rr * rr;
  /* Q(r) = -1/4 + r/5 - r^2/6 + r^3/7 - r^4/8 + r^5/9 - r^6/10, Estrin with FMA */
  double q01 = FMA(rr, 0x1.999999999999ap-3, -0.25);
  double q23 = FMA(rr, 0x1.2492492492492p-3, -0x1.5555555555555p-3);
  double q45 = FMA(rr, 0x1.c71c71c71c71cp-4, -0.125);
  double D = (r2 * r2) * FMA(r2, FMA(r2, FMA(r2, -0x1.999999999999ap-4, q45), q23), q01);
  /* sum: k ln2 + L + rh - sh/2 + ch, error terms collected in low */
  double a1 = kd * LOG_LN2H;                    /* exact */
  double s1 = a1 + Lh, e1 = Lh - (s1 - a1);     /* Fast2Sum: k = 0 or |k ln2| > |L| */
  double s2 = s1 + rh, b2 = s2 - s1, e2 = (s1 - (s2 - b2)) + (rh - b2);  /* TwoSum */
  double hs = -0.5 * sh;
  double s3 = s2 + hs, e3 = hs - (s3 - s2);     /* Fast2Sum */
  double s4 = s3 + ch, e4 = ch - (s4 - s3);     /* Fast2Sum */
  double low = (kd * LOG_LN2L + Ll) + ((rl + e1) + (e2 + e3)) + (e4 + (-0.5 * sl - rh * rl))
               + ((cl + sh * rl) + D);
  double u = s4 + low, t = low - (u - s4);      /* s4 + low = u + t exactly */
  DBG(s4, low);
  if (!(fabs_(t) > EPS * fabs_(u))) FB_LOG(x);
  if ((t > 0) != (u > 0)) u = asd(asu(u) - 1);  /* toward zero: decrement magnitude */
  return u;
}

/* cache + evaluation, compiled for FMA; only called when enabled (which implies an FMA CPU) */
static TGT double exp_fast(double x) {
  uint64_t b = asu(x); double r;
  if (cget(&c_exp, b, &r)) return r;
  return cput(&c_exp, b, exp_core(x));
}
static TGT double log_fast(double x) {
  uint64_t b = asu(x); double r;
  if (cget(&c_log, b, &r)) return r;
  return cput(&c_log, b, log_core(x));
}

double NAME_EXP(double x) { return ENABLED ? exp_fast(x) : REAL_EXP(x); }
double NAME_LOG(double x) { return ENABLED ? log_fast(x) : REAL_LOG(x); }
