# check that the dependencies are installed
pip install twine wheel


# first, let's build the project itself at all
pip install .
# or, if we are re-building it
pip install --upgrade .

# from that pont we should be able to perform a call to import from raw python
# and use the command line itself as well.


# Now, let's build the project itself
python setup.py check
# and create a distribution for it
python setup.py sdist
# as well as a wheel (weeeee)
python setup.py bdist_wheel

# now its time to upload to testpypi
twine upload --repository-url https://test.pypi.org/legacy/ dist/bioflow-<current_version>.tar.gz

# sign it
gpg --detach-sign -a dist/pyexample-0.1.0.tar.gz
# and upload the signature to the pytest:
twine upload --repository-url https://test.pypi.org/legacy/ dist/bioflow-<current_version>.tar.gz bioflow-<current_version>.tar.gz.asc

# test that pip install things properly from the testpypi:
pip install --index-url https://test.pypi.org/simple/ bioflow

# now let's send things to the actual pypi
python setup.py sdist bdist_wheel
twine check dist/bioflow-<current_version>.tar.gz
twine upload -u $PYPIUSER -p $PYPIPWD dist/bioflow-<current_version>.tar.gz bioflow-<current_version>.tar.gz.asc

# finally, let's check if it all works well
pip install bioflow