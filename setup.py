import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()


setuptools.setup(
     name='autoENRICH',
     version='0.6',
     scripts=['autoENRICH.py', 'IMPRESSION.py', 'aE_utils.py'] ,
     author="Will Gerrard",
     author_email="will.gerrard@bristol.ac.uk",
     description="Computational NMR Library",
     long_description=long_description,
   long_description_content_type="text/markdown",
     url="https://github.com/wg12385/auto-ENRICH",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)",
         "Operating System :: OS Independent",
     ],
 )
