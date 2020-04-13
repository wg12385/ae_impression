# Copyright 2020 Will Gerrard
#This file is part of autoenrich.

#autoenrich is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#autoenrich is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with autoenrich.  If not, see <https://www.gnu.org/licenses/>.

def get_unique_part(files):

	parts = []

	for file in files:
		parts.append(file.split('/')[-1].split('.')[0].split('_'))

	label_part = 0
	success = False
	while not success:
		ids = []
		for file in files:
			try:
				ids.append(file.split('/')[-1].split('.')[0].split('_')[label_part])
			except:
				label_part = -1
				success = True

		if not success:
			if len(ids) == len(set(ids)):
				success = True
			else:
				label_part += 1

	return label_part
