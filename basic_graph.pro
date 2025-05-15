QT += widgets

# Add Qt include paths
INCLUDEPATH += $$[QT_INSTALL_HEADERS]
INCLUDEPATH += $$[QT_INSTALL_HEADERS]/QtWidgets
INCLUDEPATH += $$[QT_INSTALL_HEADERS]/QtCore
INCLUDEPATH += $$[QT_INSTALL_HEADERS]/QtGui

# Add debug configuration
CONFIG += debug_and_release
CONFIG += c++20

HEADERS       = window.h \
                interpolation.h \
                approximation.h \
                graphics.h \
                tests.h

SOURCES       = main.cpp \
                window.cpp \
                interpolation.cpp \
                approximation.cpp \
                graphics.cpp

# QMAKE_CXXFLAGS += -O3 -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=2 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual
QMAKE_CXXFLAGS += -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
