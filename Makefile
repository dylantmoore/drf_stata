# Makefile for drf_plugin (C++ wrapper around drf/grf library)
#
# Targets:
#   all            - Build for local platform (darwin-arm64)
#   darwin         - Build for macOS Apple Silicon (arm64)
#   windows        - Cross-compile for Windows x86_64 (requires mingw-w64)
#   linux          - Cross-compile for Linux x86_64 (requires Linux or cross-compiler)
#   all-platforms  - Build all three targets
#   clean          - Remove all built plugins

VENDOR_DRF = vendor/drf
DRF_CORE   = $(VENDOR_DRF)/core/src
DRF_3RDPTY = $(VENDOR_DRF)/core/third_party

INCLUDES = -I$(DRF_CORE) -I$(DRF_3RDPTY) -I. -Wno-deprecated-declarations

# All drf/grf C++ source files
DRF_SRCS = \
    $(DRF_CORE)/commons/Data.cpp \
    $(DRF_CORE)/commons/DefaultData.cpp \
    $(DRF_CORE)/commons/SparseData.cpp \
    $(DRF_CORE)/commons/utility.cpp \
    $(DRF_CORE)/forest/Forest.cpp \
    $(DRF_CORE)/forest/ForestOptions.cpp \
    $(DRF_CORE)/forest/ForestPredictor.cpp \
    $(DRF_CORE)/forest/ForestPredictors.cpp \
    $(DRF_CORE)/forest/ForestTrainer.cpp \
    $(DRF_CORE)/forest/ForestTrainers.cpp \
    $(DRF_CORE)/prediction/Prediction.cpp \
    $(DRF_CORE)/prediction/PredictionValues.cpp \
    $(DRF_CORE)/prediction/RegressionPredictionStrategy.cpp \
    $(DRF_CORE)/prediction/collector/OptimizedPredictionCollector.cpp \
    $(DRF_CORE)/prediction/collector/SampleWeightComputer.cpp \
    $(DRF_CORE)/prediction/collector/TreeTraverser.cpp \
    $(DRF_CORE)/relabeling/NoopRelabelingStrategy.cpp \
    $(DRF_CORE)/sampling/RandomSampler.cpp \
    $(DRF_CORE)/sampling/SamplingOptions.cpp \
    $(DRF_CORE)/splitting/RegressionSplittingRule.cpp \
    $(DRF_CORE)/splitting/FourierSplittingRule.cpp \
    $(DRF_CORE)/splitting/factory/RegressionSplittingRuleFactory.cpp \
    $(DRF_CORE)/splitting/factory/FourierSplittingRuleFactory.cpp \
    $(DRF_CORE)/tree/Tree.cpp \
    $(DRF_CORE)/tree/TreeOptions.cpp \
    $(DRF_CORE)/tree/TreeTrainer.cpp \
    $(DRF_CORE)/analysis/SplitFrequencyComputer.cpp

# Plugin source
PLUGIN_SRC = drf_plugin_cpp.cpp

# stplugin.c — compiled as C separately per platform (different -D flags)
STPLUGIN_SRC = stplugin.c

# Output filenames
TARGET_DARWIN  = drf_plugin.darwin-arm64.plugin
TARGET_WINDOWS = drf_plugin.windows-x86_64.plugin
TARGET_LINUX   = drf_plugin.linux-x86_64.plugin

# ── Darwin (macOS arm64) ─────────────────────────────────────────────
DARWIN_CXX    = g++
DARWIN_CC     = gcc
DARWIN_CXXFLAGS = -std=c++14 -O3 -fPIC -DSYSTEM=APPLEMAC -arch arm64 $(INCLUDES)
DARWIN_CFLAGS   = -O3 -fPIC -DSYSTEM=APPLEMAC -arch arm64 -I.
DARWIN_LDFLAGS  = -bundle -lpthread

# ── Windows (x86_64, cross-compiled with mingw-w64) ─────────────────
WIN_CXX    = x86_64-w64-mingw32-g++
WIN_CC     = x86_64-w64-mingw32-gcc
WIN_CXXFLAGS = -std=c++14 -O3 -DSYSTEM=STWIN32 $(INCLUDES)
WIN_CFLAGS   = -O3 -DSYSTEM=STWIN32 -I.
WIN_LDFLAGS  = -shared -static-libstdc++ -static-libgcc -lpthread

# ── Linux (x86_64) ──────────────────────────────────────────────────
LINUX_CXX    = g++
LINUX_CC     = gcc
LINUX_CXXFLAGS = -std=c++14 -O3 -fPIC -DSYSTEM=OPUNIX $(INCLUDES)
LINUX_CFLAGS   = -O3 -fPIC -DSYSTEM=OPUNIX -I.
LINUX_LDFLAGS  = -shared -static-libstdc++ -static-libgcc -lpthread

# ── Phony targets ───────────────────────────────────────────────────
.PHONY: all darwin windows linux all-platforms clean

# Default: build for the local platform only
all: darwin

darwin: $(TARGET_DARWIN)
windows: $(TARGET_WINDOWS)
linux: $(TARGET_LINUX)
all-platforms: darwin windows linux

# ── Build rules ─────────────────────────────────────────────────────

# Darwin — compile stplugin.c as C, then link everything as C++
$(TARGET_DARWIN): $(PLUGIN_SRC) $(DRF_SRCS) $(STPLUGIN_SRC)
	$(DARWIN_CC) $(DARWIN_CFLAGS) -c $(STPLUGIN_SRC) -o stplugin.darwin.o
	$(DARWIN_CXX) $(DARWIN_CXXFLAGS) $(DARWIN_LDFLAGS) -o $@ $(PLUGIN_SRC) $(DRF_SRCS) stplugin.darwin.o
	rm -f stplugin.darwin.o

# Windows — compile stplugin.c as C with mingw, then link everything as C++
$(TARGET_WINDOWS): $(PLUGIN_SRC) $(DRF_SRCS) $(STPLUGIN_SRC)
	$(WIN_CC) $(WIN_CFLAGS) -c $(STPLUGIN_SRC) -o stplugin.windows.o
	$(WIN_CXX) $(WIN_CXXFLAGS) $(WIN_LDFLAGS) -o $@ $(PLUGIN_SRC) $(DRF_SRCS) stplugin.windows.o
	rm -f stplugin.windows.o

# Linux — compile stplugin.c as C, then link everything as C++
$(TARGET_LINUX): $(PLUGIN_SRC) $(DRF_SRCS) $(STPLUGIN_SRC)
	$(LINUX_CC) $(LINUX_CFLAGS) -c $(STPLUGIN_SRC) -o stplugin.linux.o
	$(LINUX_CXX) $(LINUX_CXXFLAGS) $(LINUX_LDFLAGS) -o $@ $(PLUGIN_SRC) $(DRF_SRCS) stplugin.linux.o
	rm -f stplugin.linux.o

clean:
	rm -f $(TARGET_DARWIN) $(TARGET_WINDOWS) $(TARGET_LINUX)
	rm -f stplugin.darwin.o stplugin.windows.o stplugin.linux.o
