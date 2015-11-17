import unittest
from os.path import join, abspath
from BioFlow.Utils.GeneralUtils.InternetIO import url_to_local_path, url_to_local_pZip, url_to_local_pGz, url_to_local


class InternetIOTest(unittest.TestCase):

    def url_to_local_path(self):
        pull_url = "http://3.bp.blogspot.com/-0N45KxrABQU/VOegbIex3BI/AAAAAAAAPzg/imoUFV-_fu0/s1600/BestPicture.jpg"
        pth = join(abspath("../../testdir/"), 'testfile.jpg')
        url_to_local_path(pull_url, pth)
        self.assertEqual(open())
