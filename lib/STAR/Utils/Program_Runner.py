
import subprocess


class Program_Runner:

    def __init__(self, cmd, scratch_dir):
        self.scratch_dir = scratch_dir
        self.executableName = cmd

    def run(self, command, cwd_dir=None):
        cmmd = command

        if not cwd_dir:
            cwd_dir = self.scratch_dir

        # print('\nRunning: ' + ' '.join(cmmd))
        p = subprocess.Popen(cmmd, cwd=cwd_dir, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, close_fds=True)
        exitCode = p.wait()

        if (exitCode == 0):
            print('\n' + ' '.join(cmmd) + ' was executed successfully, exit code was: ' +
                  str(exitCode))
        else:
            star_msg = ''
            for line in p.stdout:
                line = line.rstrip() + '\n'
                star_msg += line
            print('Error running command: ' + ' '.join(cmmd) + 'Exit Code: ' +
                  str(exitCode) + '\n\n******STAR run report******\n' + star_msg)
        return exitCode

