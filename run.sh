for nodes in {1..16}
do
  for size in {500..5000..500}
  do
    ./many.sh ${nodes} ${size} > tmp.batch
    sbatch tmp.batch
  done
done
