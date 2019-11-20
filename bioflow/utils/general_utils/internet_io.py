"""
Module responsible for retrieval of files stored on the internet.
Requires some refactoring before can be considered a library in its onw right.
"""
from os.path import abspath, join, isdir
import os
import requests
import shutil
import zipfile
import gzip
# import tarfile
import io
import requests_ftp
import hashlib


# TODO: refactor for better separation of:
#   - path correction
#   - http v.s. ftp pull selection
#   - decompression algorithm selection


def url_to_local_path(url, path, rename=None):
    """
    Copies a file from an http url to a local destination provided in path.
    Performs file-to-folder converstion
    :param url:
    :param path:
    :return:
    """
    if isdir(path) and '.zip' not in url and '.tar' not in url:
        new_path = join(path, url.split('/')[-1])

        if rename is not None:
            new_path = join(path, rename)

    if not url[:3] == 'ftp':

        r = requests.get(url, stream=True)

    else:
        # print 'debug firing'
        requests_ftp.monkeypatch_session()
        s = requests.Session()
        r = s.get(url)
        # print r.status_code
        # print r.content

    if r.status_code in ['226', 200, 226, '200']:
        if not url[:3] == 'ftp':
            with open(new_path, 'wb') as f:
                r.raw.decode_content = True
                shutil.copyfileobj(r.raw, f)
        else:
            with open(new_path, 'wb') as f:
                f.write(r.content)

    else:
        print(r.status_code)
        raise Exception(
            "Something is wrong with the url provided: %s.\n Please attempt downloading files manually" %
            url)


def url_to_local_p_zip(url, path):
    """
    Copies a file from an http url to a local folder provided  in path
    :param url:
    :param path:
    :return:
    """
    if not isdir(path):
        raise Exception("path provided %s is not a directory" % path)

    r = requests.get(url, stream=True)
    if r.status_code == 200:
        r.raw.decode_content = True
        z = zipfile.ZipFile(io.StringIO(r.content))
        z.extractall(path)  # This is unsafe as hell

    else:
        raise Exception(
            "Something is wrong with the url provided: %s.\n Please attempt downloading files manually" %
            url)


def url_to_local_p_gz(url, path):
    """
    Copies a file from an http or ftp url to a local destination provided in path
    :param url:
    :param path:
    :return:
    """
    if url[:3] == 'ftp':
        requests_ftp.monkeypatch_session()
        s = requests.Session()
        r = s.retr(url)
    else:
        r = requests.get(url, stream=True)
    if r.status_code in ['226', 200, 226, '200']:
        r.raw.decode_content = True
        f_out = open(path, 'wb')
        f_in = gzip.GzipFile(fileobj=io.StringIO((r.content)))
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
    else:
        raise Exception(
            "Something is wrong with the url provided: %s.\n Please attempt downloading files manually" %
            url)


def url_to_local(url, path, rename=None):
    """
    Copies a file from an http or ftp url to a local destination provided in path while choosing a good decompression algorithm
    so far, only gunzipped ftp url downloads and path autocompletion only for non-compressed files are supported

    :param url:
    :param path:
    :return:
    """
    if url[-2:] == 'gz':
        if rename is not None:
            raise Exception('rename unsupported for gunzipped files')
        url_to_local_p_gz(url, path)
    elif url[-3:] == 'zip':
        if rename is not None:
            raise Exception('rename unsupported for zipped files')
        url_to_local_p_zip(url, path)
    else:
        url_to_local_path(url, path, rename)


def marbach_post_proc(local_directory):
    """

    :return:
    """
    relevant_path = join(local_directory, 'Network_compendium/Tissue-specific_regulatory_networks_FANTOM5-v1/32_high-level_networks')
    for fle in os.listdir(relevant_path):
        # print fle
        if fle[-3:] == '.gz':
            f_in = join(relevant_path, fle)
            f_out = join(local_directory, fle[:-3])
            with gzip.open(f_in, 'rb') as _in, open(f_out, 'wb') as _out:
                # print _in, _out
                shutil.copyfileobj(_in, _out)

    shutil.rmtree(join(local_directory, 'Network_compendium'))


def check_hash(file_path, expected_hash, hasher):
    """
    Checks a expected_hash of a file
    :param file_path:
    :param expected_hash:
    :param hash_type:
    :return:
    """
    a_file = open(file_path, 'rb')
    block_size = 65536
    buf = a_file.read(block_size)
    while len(buf) > 0:
        hasher.update(buf)
        buf = a_file.read(block_size)
    hex_digest = hasher.hexdigest()
    return hex_digest == expected_hash


if __name__ == "__main__":
    pull_url = "http://3.bp.blogspot.com/-0N45KxrABQU/VOegbIex3BI/AAAAAAAAPzg/imoUFV-_fu0/s1600/BestPicture.jpg"
    pth = join(abspath("../../testdir/"), 'testfile.jpg')
    url_to_local_path(pull_url, pth)

    pull_url = r'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz'
    pth = abspath('../../../dumps/uniprot_sprot.dat')
    url_to_local(pull_url, pth)

    print(check_hash(pth, '439f9bf72102af7184d631d0476997d3', hashlib.md5()))
