/* portability.h */
#ifndef ELLE_PORTABILITY_H
#define ELLE_PORTABILITY_H

/* Mark an intentionally unused parameter/variable (C and C++). */
#if defined(__cplusplus)
#define ELLE_MAYBE_UNUSED [[maybe_unused]]
#else
#if defined(__GNUC__) || defined(__clang__)
#define ELLE_MAYBE_UNUSED __attribute__((unused))
#else
#define ELLE_MAYBE_UNUSED
#endif
#endif

/* One-shot suppressor for locals/results (C/C++). */
#define ELLE_UNUSED(x) (void)(x)

#endif /* ELLE_PORTABILITY_H */
