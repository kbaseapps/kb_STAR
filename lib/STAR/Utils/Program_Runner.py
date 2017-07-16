import os.path
import subprocess


class Program_Runner:

    def __init__(self, cmd, scratch_dir):
        self.scratch_dir = scratch_dir
        self.executableName = cmd

    def run(self, command, cwd=None):
        ''' options is an array of command-line parameters passed to the RQCFilter App '''
        cmmd = ' '.join(command)
        cmmd = self.executableName + ' ' + cmmd

        if not cwd:
          cwd = self.scratch_dir

        print('Running: ' + cmmd)
        p = subprocess.Popen(cmmd, cwd=cwd, shell=False)
        exitCode = p.wait()

        if (exitCode == 0):
            print(cmmd + ' was executed successfully, exit code was: ' + str(exitCode))
        else:
            raise ValueError('Error running command: ' + cmmd + '\n' +
                             'Exit Code: ' + str(exitCode))

        return exitCode

