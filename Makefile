# simple setup for a virtual environment

setup:
	@virtualenv --no-site-packages -p /usr/local/bin/python3.5 venv 	
	@venv/bin/pip install -r requirements.txt
	@echo "finished setting up the virtual environment"
	@echo "run 'source /venv/bin/activate' to start the using the virtual environment"
