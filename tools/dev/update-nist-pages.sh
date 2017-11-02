# update nist-pages
cp -r buildtemplate build
cd build
cmake -DUSE_SPHINX=ON .
make html
version=`git describe`
#branch=`git branch | grep \* | cut -d ' ' -f2`
git checkout nist-pages
cp -r html/* ../
git status
echo "Press [Enter] to commit with message: "$version"?"
read -s
git commit -a -m "$version"
#git checkout $branch
