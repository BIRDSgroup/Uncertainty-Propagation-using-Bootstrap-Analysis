#conda activate mady
B=1000

#tissues=("Whole_Blood" "Muscle_Skeletal" "Lung" "Skin_Sun_Exposed_Lower_leg" "Thyroid" "Pancreas" "Brain_Cortex" "Pituitary" "Brain_Cerebellum" "Stomach" "Brain_Caudate_basal_ganglia" "Kidney_Cortex" "Brain_Substantia_nigra" "Uterus" "Vagina" "Brain_Amygdala" "Whole_Blood_73" "Whole_Blood_237" "Muscle_Skeletal_237" "Muscle_Skeletal_73" "Lung_237" "Lung_73")

#tissues=("Whole_Blood_237" "Muscle_Skeletal_237" "Whole_Blood_73" "Muscle_Skeletal_73")

tissues=("Muscle_Skeletal_gtex" "Muscle_Skeletal_recount")

for tissue in "${tissues[@]}"
do
	echo "$tissue"
	time ./main.sh $tissue $B >> 'output.txt' 2>> 'error.txt'
	mv 'output.txt' 'error.txt' $tissue
done
#conda deactivate
