#!/usr/bin/env python

import sys
import re
import string
import os

def doit(files):


    # regular expression for ' use modulename, only: stuff, other stuff'
    # see (txt2re.com)
    use_re = re.compile("( )(use)( )((?:[a-z_][a-z_0-9]+))", 
                        re.IGNORECASE|re.DOTALL)

    module_re = re.compile("( )(module)( )((?:[a-z][a-z_0-9]+))",
                           re.IGNORECASE|re.DOTALL)

    module_proc_re = re.compile("( )(module)( )(procedure)( )((?:[a-z][a-z_0-9]+))",
                                re.IGNORECASE|re.DOTALL)

    # first parse the files and find all the module statements.  Keep a
    # dictionary of 'module name':filename.
    modulefiles = {}

    for file in files:

        f = open(file, "r")
        
        line = f.readline()

        while (line):

            # strip off the comments
            idx = string.find(line, "!")
            line = line[:idx]

            rebreak = module_re.search(line)
            rebreak2 = module_proc_re.search(line)
            if (rebreak and not rebreak2):
                modulefiles[rebreak.group(4)] = file

            line = f.readline()


        f.close()


    # go back through the files now and look for the use statements.
    # Assume only one use statement per line.  Ignore any only clauses.
    # Build a list of dependencies for the current file and output it.
    print " "
    print "# modules dependencies:"
    for file in files:

        f = open(file, "r")

        line = f.readline()
        while (line):

            # strip off the comments
            idx = string.find(line, "!")
            line = line[:idx]

            rebreak = use_re.search(line)
            if (rebreak):
                print os.path.basename(file).replace(".f90", ".o"), ':', \
                    modulefiles[rebreak.group(4)].replace(".f90", ".o")

            line = f.readline()


        f.close()
        print " "

if __name__ == "__main__":

    files = sys.argv[1:]

    doit(files)



