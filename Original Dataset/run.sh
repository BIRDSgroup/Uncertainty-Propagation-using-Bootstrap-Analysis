tissues=("Whole_Blood" "Muscle_Skeletal" "Lung" "Skin_Sun_Exposed_Lower_leg" "Thyroid" "Pancreas" "Brain_Cortex" "Pituitary" "Brain_Cerebellum" "Stomach" "Brain_Caudate_basal_ganglia" "Kidney_Cortex" "Brain_Substantia_nigra" "Uterus" "Vagina" "Brain_Amygdala")
cores=5

for tissue in "${tissues[@]}";
do
	echo "$tissue"
	time Rscript preprocessing.R $tissue $cores >> 'output.txt' 2>> 'error.txt'
done

echo "Starting for reduced samples"
tissues=("Whole_Blood" "Muscle_Skeletal" "Lung" "Skin_Sun_Exposed_Lower_leg" "Thyroid")
for tissue in "${tissues[@]}";
do
	echo "$tissue"
	time Rscript preprocessing_reduced.R $tissue >> 'output.txt' 2>> 'error.txt'
done