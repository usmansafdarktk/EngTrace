import json

file_path = 'validation_samples.jsonl'

try:
    with open(file_path, 'r', encoding='utf-8') as f:
        print(f"{'ID':<20} | {'DOMAIN':<20} | {'LEVEL':<10}")
        print("-" * 55)
        
        for line in f:
            if line.strip():  # Skip empty lines
                data = json.loads(line)
                
                # Extract fields
                sample_id = data.get('id')
                domain = data.get('domain')
                level = data.get('level')
                
                print(f"{str(sample_id):<20} | {str(domain):<20} | {str(level):<10}")

except FileNotFoundError:
    print(f"Error: The file '{file_path}' was not found.")
except Exception as e:
    print(f"An error occurred: {e}")