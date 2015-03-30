.PHONY: all help
help:

GSL_VERSION := 1.16
GSL_INSTALL_DIR := $(shell pwd)/gsl_$(GSL_VERSION)_build
GSL_TGZ := gsl-$(GSL_VERSION).tar.gz
$(GSL_TGZ):
	wget "http://gnu.mirrorcatalogs.com/gsl/$@" -O $@

.GSL_INSTALLED: $(GSL_TGZ)
	tar xzf "$<"
	mkdir -p $(GSL_INSTALL_DIR)
	cd gsl-$(GSL_VERSION) && ./configure --prefix=$(GSL_INSTALL_DIR) && make && make install
	touch "$@"

all: .GSL_INSTALLED maxent_mwe

BUILD_WHERE = build
C_SRC_ = $(shell find src -type f -name '*.c')
C_OBJ_ = $(patsubst %.c,%.o,$(subst src/,$(BUILD_WHERE)/,$(C_SRC_)))
STANDARD = -std=c99
CFLAGS := -O2 -g
$(BUILD_WHERE)/%.o: src/%.c
	@if [ ! -d $(BUILD_WHERE) ]; then mkdir -p $(BUILD_WHERE); fi
	$(CC) $(STANDARD) $(CFLAGS) -L$(GSL_INSTALL_DIR)/lib -I$(GSL_INSTALL_DIR)/include -o $@ -c $<

LIBS = -lgsl -lgslcblas -lm 

maxent_mwe: .GSL_INSTALLED $(C_OBJ_)
	$(CC) $(STANDARD) $(CFLAGS) -L$(GSL_INSTALL_DIR)/lib -o $@ $(C_OBJ_) $(LIBS)

.PHONY: clean
clean:
	rm -rf maxent_mwe $(BUILD_WHERE)
	find . -type f -name '*~' -delete

.PHONY: gsl-clean
gsl-clean:
	rm -rf .GSL_INSTALLED $(GSL_INSTALL_DIR) $(GSL_TGZ)

