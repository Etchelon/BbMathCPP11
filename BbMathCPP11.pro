TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS = -std=c++11

SOURCES += main.cpp \
    BbVector.cpp \
    BbStopWatch.cpp \
    BbMatrix.cpp

HEADERS += BbVector.hpp \
    BbStopWatch.hpp \
    BbMatrix.hpp \
    BbFactorizedMatrix.hpp \
    ColumnIterator.hpp

