Installation
============

To work, this suite of programs needs a proper installation of CMake
(http://www.cmake.org).

- on Windows it is built using Visual Studio and depends on Intel MKL (http://software.intel.com/en-us/intel-mkl/).

- on Mac OS X, Intel MKL is optional.

- on Linux it requires a development version of ATLAS (http://math-atlas.sourceforge.net/).

To compile the software do :

    Linux & Mac OS X
    ----------------
            >> cd build_directory

            build_directory is the place where you want to build OpenMEEG.
            It should be different from the OpenMEEG source directory

            >> ccmake openmeeg_source_directory

            Then choose the compilation options (ATLAS or MKL, Debug or Release ...)

            >> make

            And in order to install it :

            >> make install

    Windows
    -------
            Open the CMake.exe and configure the project
            (windows equivalent of ccmake)

            Generate the Visual Studio Solution and build it.

Quick Start
===========

Sample demo files are available at:

https://gforge.inria.fr/frs/?group_id=435

It contains scripts to compute forward solution a.k.a. lead fields
from Windows (.bat file) and on Linux or Mac OS X (bash scripts).

If OpenMEEG Python is available it can also be built from python.

Read the tutorial on line for more info on the underlying machinery:

http://openmeeg.gforge.inria.fr


Contribution rules
==================

For developpers wishing to contribute in OpenMEEG, there are several coding
rules we have fixed to have an homogeneous source code. Rules are the
following :

 - namings :
   - english is used for every class names / variable names, function names...
   - should it be in class names, functions, variables or whatever,
     abbreviations are not tolerated. All the names have to be as explicit as
     possible (the code is harder to write if you don't take advantage of
     automatic completion, but it is times easier to read for people who did
     not write this code !)
   - namespace names are lower case
   - class names start with an upper case caracter
   - variables start with a lower case caracter
   - member variables start with m_
   - global variables are prohibited
     (please contact the dev coordinator if this is necessary in your case)
   - function names do not have prefix, start with a lower case caracter and
     use upper case caracter for each new word (CamelCase : http://en.wikipedia.org/wiki/CamelCase)
   - for and if / else blocks always have curly brackets,
     even if only one call is to be done
 - english is used for the documentation of the code
 - code is documented with doxygen (http://www.stack.nl/~dimitri/doxygen/)
 - implementation is documented for complex things with or without doxygen
 - redundant include should be avoided thanks to #ifndef #define #endif in the
   begining/end of the header files.
        Ex : MyFile.h should be protected using
            #ifndef MYFILE_H
            #define MYFILE_H
            /* Code */
            #endif
 - 'using namespace' directives shall not be used in header files
 - tabs ('\t') are prohibited and are replaced with 4 spaces

This is a sample of OpenMEEG-compliant code illustrating these rules :

/**
  * \file MyClass.h
  * \author me
  * \date today
  * \brief a sample file containing MyClass class
  *
  * This sample file contains the definition of the
  * MyClass class blah blah more details...
  */
namespace myclassnamespace // use namespace if necessary
{
        /**
         * \class MyClass
         * \author me
         * \date today
         * \brief short-blah
         *
         * Detailed blah blah
          */
        class MyClass
        {
        public:
                /**
                  * \function sampleFunction
                  * \param inputValue[in] : blah
                  * \param outputValue[out] : blah
                  * \return blah blah
                  * \brief short-blah
                  *
                  * Detailed blah blah
                  */
                bool sampleFunction( int inputValue, int& outputValue )
                {
                        int counter = 0;
                        for (int i = 0; i < inputValue; ++i) {
                            counter = counter + i;
                        }
                        if(counter > 50) {
                            outputValue = counter;
                            return true;
                        }
                        return false;
                }

                /**
                  * \function content
                  * \return blah blah
                  * \brief short-blah
                  *
                  * Detailed blah blah
                  */
                int content() const { return m_content; };

        private:

            int m_content;

        };
};


