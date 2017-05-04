from cabsDock import pdb, atom
from trajectory import Trajectory
import numpy
class Trajectory_dummy(object):
	"""docstring for Trajectory_dummy"""
	def __init__(self):
		super(Trajectory_dummy, self).__init__()
		self.models = numpy.array([ ["r1", "r2", "r3"], ["r11", "r12"] ])
	def __iter__(self):
		return iter(self.models)
	def __getitem__(self, i):
		return self.models[i]
		

class Utilities(object):
	"""docstring for Utilities"""
	def __init__(self):
		super(Utilities, self).__init__()

	@staticmethod
	def perform_sequence_alignment(protein1, protein2):
		pass

	@staticmethod
	def fit_models(mobile = None, target = None, alignment = None):
		# alignment expected as [ "--ADADASDA--ADS", "ASDASDDA---DSDA"]
		# mobile and target expected as Trajectory instances
		if mobile is None or target is None:
			raise Exception("Missing input structures.")

		# get CAs from the structures.
		mobile.atoms = mobile.atoms.select("name CA")
		target.atoms = target.atoms.select("name CA")

		#extract receptors:
		chain_selection = "chain "+mobile.receptor_id[0]
		print(chain_selection)
		if len(mobile.receptor_id)>1:
			for chain_id in mobile.receptor_id:
				chain_selection += " or chain " + chain_id

		target_receptor = Trajectory(target.atoms.select(chain_selection), target.receptor_id, target.ligand_id, target.num, target.headers)
		mobile_receptor = Trajectory(mobile.atoms.select(chain_selection), mobile.receptor_id, mobile.ligand_id, mobile.num, mobile.headers)
		print("mobile receptor len ", len(mobile_receptor.atoms))
		print("mobile receptor len ", len(mobile_receptor[0].atoms))
		
		target_position = target_receptor[0].atoms.cent_of_mass()

		final_models = atom.Atoms()
		if alignment is not None: 
			# get labels of aligned residues: 
			mobile_labels, target_labels = _alignment_to_labels(alignment)
			# set of aligned target atoms
			target_reordered = Utilities._reorder(target_receptor[0], target_labels)
			for mobile_model in mobile_receptor:
				# set of aligned mobile atoms
				mobile_model_reordered = Utilities._reorder(model_receptor, mobile_labels)
				# Kabsch rotation matrix calculation for this model with ALIGNED ATOMS ONLY
				kabsch_rotation_matrix = Utilities.kabsch(mobile_model_reordered, target_reordered)
				# rotation of this model
				mobile_model_rotated = Utilities._move_and_rotate(mobile_model, target_position, kabsch_rotation_matrix)
				# storing
				final_models.append(mobile_model_rotated)

		else:
			for mobile_model, mobile_model_receptor in zip(mobile, mobile_receptor):
				print(len(mobile_model_receptor.atoms))
				kabsch_rotation_matrix = Utilities.kabsch(mobile_model_receptor, target_receptor[0])
				# rotation of this model
				mobile_model_rotated = Utilities._move_and_rotate(mobile_model, target_position, kabsch_rotation_matrix)
				# storing
				final_models+=mobile_model_rotated
		return Trajectory(final_models,mobile.receptor_id, mobile.ligand_id, mobile.num, mobile.headers) #powinno zwracac obiekt Trajectory!


	@staticmethod
	def fit_models_TEMP(mobile, target):
		# alignment expected as [ "--ADADASDA--ADS", "ASDASDDA---DSDA"]
		# mobile and target expected as Trajectory instances

		# get CAs from the structures.
		mobile.atoms = mobile.atoms.select("name CA")
		target.atoms = target.atoms.select("name CA")

		#extract receptors:
		chain_selection = "chain "+mobile.receptor_id[0]
		print(chain_selection)
		if len(mobile.receptor_id)>1:
			for chain_id in mobile.receptor_id:
				chain_selection += " or chain " + chain_id

		target_receptor = target.atoms.select(chain_selection)
		mobile_receptor = mobile.atoms.select(chain_selection)
		print("mobile receptor len ", len(mobile_receptor.atoms))
		
		target_position = target_receptor.atoms.cent_of_mass()

		final_models = atom.Atoms()
		for mobile_model, mobile_model_receptor in zip(mobile.get_models(), mobile_receptor.get_models()):
			print(len(mobile_model_receptor.atoms))
			kabsch_rotation_matrix = Utilities.kabsch(mobile_model_receptor, target_receptor)
			# rotation of this model
			mobile_model_rotated = Utilities._move_and_rotate(mobile_model, target_position, kabsch_rotation_matrix)
			# storing
			final_models+=mobile_model_rotated
		return final_models


	@staticmethod
	def _move_and_rotate(model = None, translation_vector = None, rotation_matrix = None):
		# moves and rotates the model atoms
		if model is None or translation_vector is None or rotation_matrix is None:
			raise Exception("Missing arguments.")
		model_center_of_mass = model.atoms.cent_of_mass()
		# (1) move to the origin, (2) rotate with the rotation matrix, (3) move back to original position and move with translation vector
		return model.atoms.center_at_origin().rotate(rotation_matrix).move(model_center_of_mass + translation_vector)

	@staticmethod
	def kabsch(model1, model2, concentric=False):
		#Computes the rotation matrix for best fit between two models.
		model1_atoms = model1.atoms
		model2_atoms = model2.atoms
		if len(model1_atoms) is not len(model2_atoms):
			raise IndexError('Atom sets have different length: %i != %i' % (len(model1_atoms), len(model2_atoms)))
		t = model2_atoms.to_matrix()
		q = model1_atoms.to_matrix()
		if not concentric:
			t = numpy.subtract(t, numpy.average(t, 1))
			q = numpy.subtract(q, numpy.average(q, 1))
		v, s, w = numpy.linalg.svd(numpy.dot(t, q.T))
		d = numpy.identity(3)
		if numpy.linalg.det(numpy.dot(w.T, v.T)) < 0:
			d[2, 2] = -1
		return numpy.matrix(numpy.dot(numpy.dot(w.T, d), v.T))


	
	@staticmethod
	def _alignment_to_labels(alignment = None):
		# alignment expected as [ "--ADADASDA--ADS", "ASDASDDA---DSDA"]
		if alignment is None or len(alignment[0]) != len(alignment[1]):
			raise Exception
		else:
			# generate a lists of residue LABELS of aligned residues
			mobile_labels = []
			target_labels = []
			mobile_residue_label = 0
			target_residue_label = 0
			for mobile_residue_flag, target_residue_flag in zip(alignment[0], alignment[1]):
				mobile_residue_flag_bool = mobile_residue_flag is not "-"
				target_residue_flag_bool = target_residue_flag is not "-"
				if mobile_residue_flag_bool:
					if target_residue_flag_bool:
						mobile_labels.append(mobile_residue_label)
						target_labels.append(target_residue_label)
						target_residue_label += 1
					mobile_residue_label += 1
				else:
					if target_residue_flag_bool:
						target_residue_label += 1
		return mobile_labels, target_labels

	@staticmethod
	def _reorder(_list, order):
		return [_list[item] for item in order]

	@staticmethod
	def filter_models(trajectories = [], _filter = "default", N_from_each = 100):
		if _filter is not "default":
			raise Exception("No filter defined.")
		else:
			filtered_models = []
			filtered_atoms = atom.Atoms()
			for trajectory in trajectories:
				models_list = trajectory.get_models()
				print("models", len(models_list))
				models_list_filtered = [ model for model in models_list if model.e_int < 0] 
				print("models_list_filtered", len(models_list_filtered))
				#print( [model.e_int for model in models_list_filtered] )
				models_list_filtered_sorted = sorted(models_list_filtered, key = lambda model: model.e_int)
				print(models_list_filtered_sorted)
				if len(models_list_filtered_sorted) <= N_from_each:
					filtered_models += models_list_filtered_sorted
				else:
					filtered_models += models_list_filtered_sorted[:N_from_each]
				model_num = 1
				for model in filtered_models:
					filtered_atoms += model.atoms.set_model_number(model_num)
					model_num += 1
		print("filtered models", len(filtered_models))
		output_file = open("sa_top_1000.pdb", 'w')
		output_file.write(filtered_atoms.make_pdb())
		output_file.close()
		return filtered_models, filtered_atoms


	@staticmethod
	def rmsd(model1, model2, selection = None):
		r = 0
		for a1, a2 in zip(model1.atoms, model2.atoms):
				r += a1.dist2(a2)
		return numpy.sqrt(r/len(model1.atoms))



# alignment = ["ABC", "C-D"]
# mlabels, tlabels = Utilities._alignment_to_labels(alignment)
# t = Trajectory_dummy()
# print(mlabels)
# for residue in Utilities._reorder(t[0], mlabels):
# 	print residue

# print(t[0])
#Utilities.fit_models()
#dsa


