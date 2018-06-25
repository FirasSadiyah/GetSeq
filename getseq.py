import sys

if (sys.version_info < (3, 0)):
    print("Python 3.0 or a recent version is required.")
else:
    import getseq_cmd

    getseq_cmd.getseq()
