
# Simple get yes or no from user
def yes_or_no():
	reply = str(input(' (y/n): ')).lower().strip()
	if len(reply) == 0:
		return True
	if reply[0] == 'y':
		return True
	if reply[0] == 'n':
		return False
	else:
		return False
