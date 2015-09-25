
header-links: remove-links
	HEADERLIST="$(top_srcdir)/config.h $(top_srcdir)/base/*.h $(top_srcdir)/util/*.h"; \
	for h in $$HEADERLIST; do \
	  BASENAME=`basename $$h`; \
	  $(LN_S) $$h $$BASENAME;\
	done

remove-links: 
	rm -f *.h


all: all-am header-links

clean: clean-am remove-links
distclean: distclean-am remove-links
