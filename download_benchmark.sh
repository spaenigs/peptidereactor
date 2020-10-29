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

#HOST=172.16.103.199
#NAME=ubuntu

DATASET_NAMES=(
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
  "aip_antiinflam"
)

HOST=172.16.103.203
NAME=ubuntu

#DATASET_NAMES=(
#  "cpp_mixed"
#  "amp_gonzales"
#  "cpp_sanders"
#  "hiv_bevirimat"
#  "tce_zhao"
#  "amp_fernandes"
#  "hiv_nfv"
#  "hiv_v3"
#)

#HOST=137.248.121.201
#NAME=spaenigs


for i in "${DATASET_NAMES[@]}"
do
	echo $i
	scp -r $NAME@$HOST:/home/$NAME/peptidereactor/data/$i/benchmark \
	  /media/spaenigs/565f856e-2d16-4784-a91c-c36edf56faf4/peptidereactor/data/$i/
	scp -r $NAME@$HOST:/home/$NAME/peptidereactor/data/$i/vis \
	  /media/spaenigs/565f856e-2d16-4784-a91c-c36edf56faf4/peptidereactor/data/$i/
done