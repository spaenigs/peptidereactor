#DATASET_NAMES=(
#  "amp_csamp"
#  "acp_iacp"
#  "cpp_mlcppue"
#  "acp_anticp"
#  "atb_antitbp"
#  "hiv_lpv"
#  "acp_mlacp"
#  "hiv_abc"
#  "hiv_azt"
#  "hiv_d4t"
#  "hiv_ddi"
#  "hiv_3tc"
#  "amp_iamp2l"
#  "hem_hemopi"
#  "aip_aippred"
#)

#HOST=172.16.103.179
#NAME=ubuntu

#DATASET_NAMES=(
#  "hiv_apv"
#  "hiv_dlv"
#  "hiv_efv"
#  "hiv_rtv"
#  "hiv_nvp"
#  "hiv_idv"
#  "hiv_sqv"
#  "amp_antibp"
#  "cpp_cppredfl"
#  "hiv_protease"
#  "cpp_kelmcpp"
#  "avp_avppred"
#  "cpp_cellppdmod"
#  "pip_pipel"
#  "hiv_sqv"
#  "atb_iantitb"
#)
#
#HOST=172.16.103.199
#NAME=ubuntu

#DATASET_NAMES=(
#  "avp_amppred"
#  "cpp_cellppd"
#  "nep_neuropipred"
#  "cpp_mlcpp"
#  "amp_antibp2"
#  "bce_ibce"
#  "amp_modlamp"
#  "afp_amppred"
#  "afp_antifp"
#  "isp_il10pred"
#  "ace_vaxinpad"
#  "aip_antiinflam"
#)
#
#HOST=172.16.103.203
#NAME=ubuntu

DATASET_NAMES=(
  "cpp_mixed"
  "amp_gonzales"
  "cpp_sanders"
  "hiv_bevirimat"
  "tce_zhao"
  "amp_fernandes"
  "hiv_nfv"
  "hiv_v3"
)

HOST=137.248.121.201
NAME=spaenigs

TARGET_DIR=/media/spaenigs/565f856e-2d16-4784-a91c-c36edf56faf4/peptidereactor/

for i in "${DATASET_NAMES[@]}"
do

	echo $i

	## 1.
	# a) delete all, sequence_based, structure_based
	rm -r $TARGET_DIR/data/$i/csv/all
	rm -r $TARGET_DIR/data/$i/csv/sequence_based
	rm -r $TARGET_DIR/data/$i/csv/structure_based
	# b) download the 3 folders
	scp -r $NAME@$HOST:/home/$NAME/peptidereactor/data/$i/csv/sequence_based/ $TARGET_DIR/data/$i/csv/
	scp -r $NAME@$HOST:/home/$NAME/peptidereactor/data/$i/csv/structure_based/ $TARGET_DIR/data/$i/csv/
	scp -r $NAME@$HOST:/home/$NAME/peptidereactor/data/$i/csv/all/ $TARGET_DIR/data/$i/csv/

	## 3.
	# a) delete misc/benchmark on server
	ssh $NAME@$HOST "rm -r /home/$NAME/peptidereactor/data/$i/misc/benchmark/"
	# b) upload misc/benchmark dir to server
	scp -r $TARGET_DIR/data/$i/misc/benchmark/ $NAME@$HOST:/home/$NAME/peptidereactor/data/$i/misc/

done