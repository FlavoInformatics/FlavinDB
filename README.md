FlavinDB
========

## Setup
1. Install python3.6 and virtualenv
    * Mac
        1. Install brew
            * `/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`
        2. Install python3.6
            * `brew update`
            * `brew install python3.6`
        3. Install virtualenv
            * `pip install virtualenv`
2. Setup the virtualenv
    * `virtualenv -p /usr/local/bin/python3.6 venv`
3. Load the virtualenv
    * `source venv/bin/activate`
3. Install the requirements
    * `pip install -r requirements.txt`