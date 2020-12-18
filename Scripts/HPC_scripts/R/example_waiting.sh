
echo "In"
for i in {1..10}; do
  sleep 20 &
  pids[${i}]=$!
done

now=$(date)
echo "$now"

# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done

now=$(date)
echo "$now"
