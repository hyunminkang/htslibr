PACKAGE_DIR=htslibr

.PHONY= install clean

install: htslibr_1.0.tar.gz

$(PACKAGE_DIR):
	Rscript -e "Rcpp::Rcpp.package.skeleton('htslibr')"

$(PACKAGE_DIR)/src/htslib/hts.c: $(PACKAGE_DIR)	
	git submodule add https://github.com/samtools/htslib.git $(PACKAGE_DIR)/src/htslib

$(PACKAGE_DIR)/src/Makevars: $(PACKAGE_DIR) Makevars
	cp Makevars $</src/

$(PACKAGE_DIR)/src/util.cpp: $(PACKAGE_DIR) util.cpp
	cp util.cpp $</src/

$(PACKAGE_DIR)/inst/include/hts.h: $(PACKAGE_DIR)
	mkdir -p $</inst/include
	cp $(PACKAGE_DIR)/src/htslib/htslib/*.h $</inst/include
	
htslibr_1.0.tar.gz: $(PACKAGE_DIR)/src/htslib/hts.c $(PACKAGE_DIR)/src/Makevars $(PACKAGE_DIR)/inst/include/hts.h $(PACKAGE_DIR)/src/util.cpp
	cd $(PACKAGE_DIR) && Rscript -e "devtools::build()"

# install: $(PACKAGE_DIR) $(PACKAGE_DIR)/src/htslib htslibr_1.0.tar.gz

clean:
	rm -rf $(PACKAGE_DIR)
