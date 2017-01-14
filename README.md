FlavinDB
========

## Setup
1. Install python3.6 and virtualenv
    * Mac
        1. Install Homebrew
            * `/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`
        2. Install python3.6
            * `brew update`
            * `brew install python3`
        3. Install virtualenv
            * `pip3 install virtualenv`
2. Setup the virtualenv
    * `virtualenv -p /usr/local/bin/python3 venv`
3. Load the virtualenv
    * `source venv/bin/activate`
3. Install the requirements
    * `pip3 install -r requirements.txt`
    * NOTE for OS X users: install `libpng` and `freetype` via brew if you get a `Command "python setup.py egg_info" failed with error code 1` error when trying to install matplotlib. See http://stackoverflow.com/questions/9829175/pip-install-matplotlib-error-with-virtualenv for more details.