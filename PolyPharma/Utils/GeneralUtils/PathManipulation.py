__author__ = 'ank'

import os

def mkdir_recursive(path):
    sub_path = os.path.dirname(path)
    print 'subpath: %s' % sub_path
    if not os.path.exists(sub_path):
        mkdir_recursive(sub_path)
    if not os.path.exists(path):
        print 'path %s does not exist yet, looks like file: %s' % (path, '.' in path.split('/')[-1][-5:]),
        if not '.' in path.split('/')[-1][-5:]:  # should be able to supress specific file creation ...
            os.mkdir(path)
            print '; created'
        else:
            print '; creation skipped'


