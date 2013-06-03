TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS = -std=c++11

SOURCES += main.cpp \
    ludecomp.cpp \
    BbVector.cpp \
    BbStopWatch.cpp \
    BbMatrix.cpp

HEADERS += \
    ludecomp.hpp \
    BbVector.hpp \
    BbStopWatch.hpp \
    BbMatrix.hpp \
    BbFactorizedMatrix.hpp \
    ColumnIterator.hpp

