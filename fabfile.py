"""
Used for two things: 1. Locally for git workflow, 2. Update remote servers.
"""
from __future__ import with_statement
from fabric.api import *
from fabric.colors import *
from fabric.context_managers import cd, lcd
import os.path

HOST = '' # e.g. seb@bla.com
REMOTE_BASE_DIR = '/webapps/seb_django/www'  # absolute path, where project/repo lives
REMOTE_ERR_FILE = '/webapps/seb_django/logs/00UPDATE_203341.err'  # absolute path
REMOTE_LOG_FILE = '/webapps/seb_django/logs/00UPDATE_203341.log'  # absolute path
REPO_NAME = 'genomics'  # basename of project
REPO_URL = 'git@github.com:sschmeier/genomics.git'  # e.g. github url
REPO_BRANCH = 'master'  # this is the branch to clone on servers


@hosts('seb@vm010865.massey.ac.nz', 'seb@vm010944.massey.ac.nz') # only for deploy
def logs():
    """ Reading Massey server log-files and print to stdout."""
    puts(yellow("[Reading log-file]"))
    run("cat %s" % REMOTE_LOG_FILE)
    puts(yellow("[Reading err-file]"))
    run("cat %s" % REMOTE_ERR_FILE)


@hosts('seb@vm010865.massey.ac.nz', 'seb@vm010944.massey.ac.nz') # only for deploy
def deploymassey(activate_env=True, conda=None):
    """ Deploy project to Massey servers."""
    remote_dir = os.path.abspath(os.path.join(REMOTE_BASE_DIR, REPO_NAME))
    
    if activate_env:
        if conda:
            puts(yellow("[Activate conda env]"))
            run('source activate %s' %(conda))
        else:
            puts(yellow("[Activate env through ~/bin/activate]"))
            run('source ~/bin/activate')

    with settings(warn_only=True):
        if run("test -d %s" % (remote_dir)).failed:
            puts(red("[Repo %s does not exist on remote at: %s]" % (REPO_NAME, remote_dir)))
            with cd(REMOTE_BASE_DIR):
                run("git clone -b %s %s %s" % (REPO_BRANCH, REPO_URL, REPO_NAME))

    puts(yellow("[Write logs]"))
    run("echo '-----------------------------' > %s" % REMOTE_ERR_FILE)
    run("echo `date` >> %s" % REMOTE_ERR_FILE)
    run("echo '-----------------------------' >> %s" % REMOTE_ERR_FILE)
    run("echo '-----------------------------' > %s" % REMOTE_LOG_FILE)
    run("echo `date` >> %s" % REMOTE_LOG_FILE)
    run("echo '-----------------------------' >> %s" % REMOTE_LOG_FILE)

    puts(yellow("[Update repo: %s]" % REPO_NAME))
    with cd(remote_dir):
        run("git pull origin %s >> %s 2>> %s" %
            (REPO_BRANCH, REMOTE_LOG_FILE, REMOTE_ERR_FILE))


def deploy(br='master', bump="patch"):
    """Deploy master to remotes.
    - make latexpdf
    - copy pdf to _static.
    - commit changes.
    - bump2version
    - push master to origin (github)

    Keyword arguments:
    br -- the branch that should be pushed
    bump -- patch, minor, major

    Usage:
    fab deploy
    """
    # co branch
    puts(yellow("[Checkout branch %s]"%(br)))
    local("git checkout %s"%(br))

    # create new pdf
    puts(yellow('[Make latexpdf]'))
    local("make latexpdf > /tmp/latex")

    puts(yellow('[Copy pdf]'))
    local("cp build/latex/*.pdf source/_static/")

    puts(yellow('[git stage/commit changes]'))
    local("git add -u")
    local('git commit -m "added new pdf"' %(msg))
    
    # push changes to gitlab
    puts(yellow("[Bump version]"))
    local("bump2version %s"%(bump))

    # push changes to github
    puts(yellow("[Push %s to github]"%(br)))
    local("git push origin %s"%(br))
