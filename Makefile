BUILD	?= /tmp/lse

BENCH	= $(BUILD)/bench

CC	= gcc
LD	= gcc
GDB	= gdb
MK	= mkdir -p
RM      = rm -rf

CFLAGS	= -std=gnu99 -Wall -O3 -flto=auto -g3 -pipe

CFLAGS  += -fno-math-errno \
	   -ffinite-math-only \
	   -fno-signed-zeros \
	   -fno-trapping-math \
	   -fno-associative-math \
	   -fno-reciprocal-math \
	   -ffp-contract=fast

LFLAGS	= -lm

OBJS	= bench.o lse.o lfg.o

BENCH_OBJS = $(addprefix $(BUILD)/, $(OBJS))

all: $(BENCH)

$(BUILD)/%.o: %.c
	@ echo "  CC    " $<
	@ $(MK) $(dir $@)
	@ $(CC) -c $(CFLAGS) -MMD -o $@ $<

$(BENCH): $(BENCH_OBJS)
	@ echo "  LD    " $(notdir $@)
	@ $(LD) $(CFLAGS) -o $@ $^ $(LFLAGS)

run: $(BENCH)
	@ echo "  RUN	" $(notdir $<)
	@ $<

debug: $(BENCH)
	@ echo "  GDB	" $(notdir $<)
	@ $(GDB) $<

clean:
	@ echo "  CLEAN "
	@ $(RM) $(BUILD)

include $(wildcard $(BUILD)/*.d)

