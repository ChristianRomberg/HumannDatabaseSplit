# assume we can align 70% of the sequences with 10% of the database
quantile = 0.7
quantile_value = 0.1

# percent of time saved shouldn't change with these, but filling proper values will yield a rough estimate of the time
# a two-staged approach would take
db_size=1
time_spent_total=1
num_samples=1

time_spent_first_stage = num_samples * (db_size * quantile_value) * time_spent_total 
time_spent_second_stage = (num_samples * (1-quantile)) * (db_size * (1-quantile_value) * time_spent_total

time_spent_both_stages = num_samples * (db_size * quantile_value) * time_spent_total + num_samples * (1-quantile) * (db_size * (1-quantile_value) * time_spent_total
	 = (num_samples * time_spent_total * db_size) * ((quantile_value) + ((1-quantile) * (1-quantile_value)))