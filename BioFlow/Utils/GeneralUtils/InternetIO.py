"""
Module responsible for retrieval of files stored on the internet.
Requires some refactoring before can be considered a library in its onw right.
"""
__author__ = 'ank'

from os.path import abspath, join, isdir
import requests
import shutil
import zipfile
import gzip
# import tarfile
import StringIO
import requests_ftp


# TODO: refactor for better separation of:
#   - path correction
#   - http v.s. ftp pull selection
#   - decompression algorithm selection


def url_to_local_path(URL, path):
    """
    Copies a file from an http URL to a local destination provided in path.
    Performs file-to-folder converstion
    :param URL:
    :param path:
    :return:
    """
    if isdir(path) and '.zip' not in URL and '.tar' not in URL:
        path = join(path, URL.split('/')[-1])
    r = requests.get(URL, stream=True)
    if r.status_code == 200:
        with open(path, 'wb') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)
    else:
        raise Exception("Something is wrong with the url provided: %s.\n Please attempt downloading files manually" % URL)


def url_to_local_pZip(URL, path):
    """
    Copies a file from an http URL to a local folder provided  in path
    :param URL:
    :param path:
    :return:
    """
    if not isdir(path):
        raise Exception("path provided %s is not a directory" % path)
    r = requests.get(URL, stream=True)
    if r.status_code == 200:
        r.raw.decode_content = True
        z = zipfile.ZipFile(StringIO.StringIO(r.content))
        z.extractall(path)  #This is unsafe as hell
    else:
        raise Exception("Something is wrong with the url provided: %s.\n Please attempt downloading files manually" % URL)


def url_to_local_pGz(URL, path):
    """
    Copies a file from an http or ftp URL to a local destination provided in path
    :param URL:
    :param path:
    :return:
    """
    if URL[:3] =='ftp':
        requests_ftp.monkeypatch_session()
        s = requests.Session()
        r = s.retr(URL)
    else:
        r = requests.get(URL, stream=True)
    if r.status_code in ['226', 200, 226, '200']:
        r.raw.decode_content = True
        f_out = open(path, 'wb')
        f_in = gzip.GzipFile(fileobj=StringIO.StringIO((r.content)))
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
    else:
        raise Exception("Something is wrong with the url provided: %s.\n Please attempt downloading files manually" % URL)


def url_to_local(URL, path):
    """
    Copies a file from an http or ftp URL to a local destination provided in path while choosing a good decompression algorithm
    so far, only gunzipped ftp URL downloads and path autocompletion only for non-compressed files are supported

    :param URL:
    :param path:
    :return:
    """
    if URL[-2:] == 'gz':
        url_to_local_pGz(URL, path)
    elif URL[-3:] == 'zip':
        url_to_local_pZip(URL, path)
    else:
        url_to_local_path(URL, path)


if __name__ == "__main__":
    pull_url = "http://3.bp.blogspot.com/-0N45KxrABQU/VOegbIex3BI/AAAAAAAAPzg/imoUFV-_fu0/s1600/BestPicture.jpg"
    pth = join(abspath("../../testdir/"), 'testfile.jpg')
    url_to_local_path(pull_url, pth)

    pull_url = r'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'
    pth = abspath('../../testdir/uniprot_sprot.dat')
    url_to_local(pull_url, pth)