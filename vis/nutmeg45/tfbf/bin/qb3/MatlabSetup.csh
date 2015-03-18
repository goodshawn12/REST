#!/bin/csh

setenv MCR_ROOT $HOME/bin/MCR/v74
setenv MATLAB_SHELL /bin/bash

switch (`arch`)
case i686
   setenv LD_LIBRARY_PATH $MCR_ROOT/runtime/glnx86:$MCR_ROOT/sys/os/glnx86:$MCR_ROOT/sys/java/jre/glnx86/jre1.5.0/lib/i386/native_threads:$MCR_ROOT/sys/java/jre/glnx86/jre1.5.0/lib/i386/client:$MCR_ROOT/sys/java/jre/glnx86/jre1.5.0/lib/i386:$MCR_ROOT/bin/glnx86
   breaksw
case x86_64
   setenv MOSTYPE glnxa64
   setenv MARCH amd64
#   setenv LD_LIBRARY_PATH $MCR_ROOT/runtime/glnxa64:$MCR_ROOT/sys/os/glnxa64:$MCR_ROOT/sys/java/jre/glnxa64/jre1.4.2/lib/amd64/native_threads:$MCR_ROOT/sys/java/jre/glnxa64/jre1.4.2/lib/amd64/client:$MCR_ROOT/sys/java/jre/glnxa64/jre1.4.2/lib/amd64:$MCR_ROOT/bin/glnxa64
   setenv LD_LIBRARY_PATH $MCR_ROOT/runtime/glnxa64:$MCR_ROOT/sys/os/glnxa64:$MCR_ROOT/sys/java/jre/glnxa64/jre1.5.0/lib/amd64/native_threads:$MCR_ROOT/sys/java/jre/glnxa64/jre1.5.0/lib/amd64/client:$MCR_ROOT/sys/java/jre/glnxa64/jre1.5.0/lib/amd64:$MCR_ROOT/bin/glnxa64
   breaksw
endsw

setenv XAPPLRESDIR $MCR_ROOT/X11/app-defaults

unsetenv MCR_ROOT
unsetenv MARCH
