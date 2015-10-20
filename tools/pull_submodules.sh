# Obtain them, if you didn't do a recursive clone
git submodule update --init --recursive

# Update them
git submodule foreach git pull origin master
git pull --recurse-submodules
