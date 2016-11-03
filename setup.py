from setuptools import setup


setup(
    name='gdemux',
    version='0.0.1',
    url='http://github.com/brwnj/gdemux',
    license='Unlicense',
    author='Joe Brown',
    author_email='joe.brown@pnnl.gov',
    description='Almost pointless demultiplexer',
    long_description=__doc__,
    py_modules=['gdemux'],
    install_requires=[
        'click',
        'pandas'
    ],
    entry_points='''
        [console_scripts]
        gdemux=gdemux:group_demux
    '''
)
