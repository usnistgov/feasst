from skbuild import setup

setup(
#    name="feasst",
#    packages=["src/feasst"],
)

# I struggled to get scikit-build-core to work with libraries, but the following setuptools guide seems to work
# from https://gist.github.com/kprussing/db21614ca5b51cedff07dfb70059f280
# https://stackoverflow.com/questions/70044257/packaging-executable-shared-library-and-python-bindings-not-finding-library
