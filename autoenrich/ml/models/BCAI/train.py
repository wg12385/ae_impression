## Copyright (c) 2017 Robert Bosch GmbH
## All rights reserved.
##
## This source code is licensed under the MIT license found in the
## LICENSE file in the root directory of this source tree.

import os
import torch
torch.manual_seed(0)
import torch.nn as nn
import torch.nn.functional as F
from torch import optim
from torch.utils.data import DataLoader
from torch import autograd
#from graph_transformer import *
import json
import time
import numpy as np
import math
import argparse
import sys
from tqdm import tqdm

np.set_printoptions(threshold=sys.maxsize)

#####################################################################
#
# Training script
#
#####################################################################

def loss(y_pred, y, x_bond):
	NUM_BOND_ORIG_TYPES = 8
	champs_loss = False

	y_pred[torch.abs(y_pred)<0.000001] = 0.0
	y_pred_pad = torch.cat([torch.zeros(y_pred.shape[0], 1, y_pred.shape[2], device=y_pred.device), y_pred], dim=1)
	# Note: The [:,:,1] below should match the num_bond_types[1]*final_dim in graph transformer
	y_pred_scaled = y_pred_pad.gather(1,x_bond[:,:,1][:,None,:])[:,0,:] * y[:,:,2] + y[:,:,1]

	abs_dy = (y_pred_scaled - y[:,:,0]).abs()
	loss_bonds = (x_bond[:,:,0] > 0)
	abs_err = abs_dy.masked_select(loss_bonds & (y[:,:,3] > 0)).sum()
	mabs_err = abs_dy.masked_select(loss_bonds & (y[:,:,3] > 0)).mean()

	type_dy = [abs_dy.masked_select(x_bond[:,:,0] == i) for i in range(1,NUM_BOND_ORIG_TYPES+1)]
	if champs_loss:
		type_err = torch.cat([t.sum().view(1) for t in type_dy], dim=0)
		type_cnt = torch.cat([torch.sum(x_bond[:,:,0] == i).view(1) for i in range(1,NUM_BOND_ORIG_TYPES+1)])
	else:
		type_err = torch.tensor([t.sum() for t in type_dy])
		type_cnt = torch.tensor([len(t) for t in type_dy])

	return abs_err, type_err, type_cnt, mabs_err


def epoch(loader, model, opt=None, lr=0.001):
	device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
	para_model = model.to(device)

	NUM_BOND_ORIG_TYPES = 8
	MAX_BOND_COUNT = 500
	log_interval = 1
	batch_chunk = 1
	batch_size = 10
	# Dont use this, i suspect not having some of every cpl type is an issue
	champs_loss = False
	warmup_step = 5000
	clip = 0.25
	APEX_AVAILABLE = False

	global train_step
	train_step = 0
	start_epoch = 0
	model.eval() if opt is None else model.train()
	dev = next(model.parameters()).device
	abs_err, type_err, type_cnt, mabs_err = 0.0, torch.zeros(NUM_BOND_ORIG_TYPES), torch.zeros(NUM_BOND_ORIG_TYPES, dtype=torch.long), 0.0
	log_interval = log_interval

	with torch.enable_grad() if opt else torch.no_grad():
		batch_id = 0
		total_loss = torch.zeros(NUM_BOND_ORIG_TYPES)
		iter = 0
		with autograd.detect_anomaly():
			with tqdm(total=len(loader), desc='tr loss: ', leave=False) as pbar:
				for x_idx, x_atom, x_atom_pos, x_bond, x_bond_dist, x_triplet, x_triplet_angle, y in loader:
					x_atom, x_atom_pos, x_bond, x_bond_dist, x_triplet, x_triplet_angle, y = \
						x_atom.to(dev), x_atom_pos.to(dev), x_bond.to(dev), x_bond_dist.to(dev), \
						x_triplet.to(dev), x_triplet_angle.to(dev), y.to(dev)

					x_bond, x_bond_dist, y = x_bond[:, :MAX_BOND_COUNT], x_bond_dist[:, :MAX_BOND_COUNT], y[:,:MAX_BOND_COUNT]
					iter += 1
					if opt:
						# Put this here so that the batch_chunk setting will work
						opt.zero_grad()

						# Perform cutout on the molecule (i.e., for a large molecule, randomly remove an atom and its nearest
						# neighbor; then remove all bonds/triplets related to this atom)
						#x_atom, x_atom_pos, x_bond, x_bond_dist, x_triplet, x_triplet_angle, _ = \
						#	subgraph_filter(x_atom, x_atom_pos, x_bond, x_bond_dist, x_triplet, x_triplet_angle, args)

					if batch_chunk > 1:
						mbsz = batch_size // batch_chunk
						b_abs_err = 0
						b_type_err = 0
						b_type_cnt = 0
						types_cnt = sum([(x_bond[:,:,0] == i).sum() for i in range(1,NUM_BOND_ORIG_TYPES+1)])
						for i in range(batch_chunk):
							mini = slice(i*mbsz,(i+1)*mbsz)
							y_pred_mb, _ = para_model(x_atom[mini], x_atom_pos[mini], x_bond[mini], x_bond_dist[mini],
													  x_triplet[mini], x_triplet_angle[mini])
							mb_abs_err, mb_type_err, mb_type_cnt = loss(y_pred_mb, y[mini], x_bond[mini])
							b_abs_err += mb_abs_err.detach()       # No need to average, as it's sum
							b_type_err += mb_type_err.detach()
							b_type_cnt += mb_type_cnt.detach()
							mb_raw_loss = mb_abs_err / types_cnt.float()
							if champs_loss:
								raise ValueError("CHAMPS loss not supported yet with batch_chunk mode")
							if APEX_AVAILABLE:
								with amp.scale_loss(mb_raw_loss, opt) as scaled_loss:
									scaled_loss.backward()
							else:
								mb_raw_loss.backward()
					else:
						y_pred, _ = para_model(x_atom, x_atom_pos, x_bond, x_bond_dist, x_triplet, x_triplet_angle)
						#print(y_pred.cpu().detach().numpy())
						b_abs_err, b_type_err, b_type_cnt, b_mabs = loss(y_pred, y, x_bond)

					abs_err += b_abs_err.detach()
					type_err += b_type_err.detach()
					type_cnt += b_type_cnt.detach()
					mabs_err += b_mabs.detach()
					batch_id += 1
					total_loss += b_type_err / b_type_cnt.float()

					if opt:
						train_step += 1
						if train_step <= warmup_step:
							# It does this
							curr_lr = lr * train_step / warmup_step
							opt.param_groups[0]['lr'] = curr_lr
							#print('\tLR: ', curr_lr)
						elif args.scheduler == 'cosine':
							# Not this
							scheduler.step(train_step)

						if batch_chunk == 1:
							# it does this
							raw_loss = b_abs_err/b_type_cnt.sum()
							raw_loss.backward()

						opt.step()

						#print('\tLOSS = ', raw_loss.detach().numpy())
						if batch_id % log_interval == 0:
							# Dont know what this is but it does it
							avg_loss = torch.log(total_loss / log_interval).mean().item()
							#logging(f"Epoch {ep:2d} | Step {train_step} | lr {opt.param_groups[0]['lr']:.7f} | Error {avg_loss:.5f}")

							total_loss = 0

					pbar.update(1)
					pbar.set_description("batch MAE: {0:<10.4f}".format(b_mabs))

	torch.cuda.empty_cache()
	return abs_err / type_cnt.sum(), torch.log(type_err / type_cnt.float()).mean(), torch.log(type_err / type_cnt.float()), mabs_err/iter






















##
