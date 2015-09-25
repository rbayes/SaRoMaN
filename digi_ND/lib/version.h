#ifndef __DIGI_VERSION_H__
#define __DIGI_VERSION_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
__BEGIN_DECLS


#define DIGI_VERSION "v0r1p0"

extern const char * digi_version;

__END_DECLS

#endif /* __DIGI_VERSION_H__ */
