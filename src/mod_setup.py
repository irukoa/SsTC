import re
from glob import glob

headerln = re.compile("!###MOD HEADERS")
procsln  = re.compile("!###MOD PROCEDURES")

#Open unmmoded source code.
with open('./src/SsTC.F90', "r") as unmodded_SsTC:
    unmoddedln = unmodded_SsTC.readlines()

#Create temporary source code with appropiate "include ..." insertions.
with open('./src/SsTC_mod.F90', "w") as modded_SsTC:
    for line in unmoddedln:
        for match in re.finditer(headerln, line):
            for filename in glob('src/**/headers.inc', recursive=True):
              print("Including header file: " + filename)
              line = line + "  include " + '"' + filename + '"\n'
        for match in re.finditer(procsln, line):
            for filename in glob('src/**/procedures.inc', recursive=True):
              print("Including procedure list file: " + filename)
              line = line + "  include " + '"' + filename + '"\n'
        modded_SsTC.write(line)
