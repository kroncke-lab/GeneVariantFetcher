import csv
from utils.variant_normalizer import VariantNormalizer

# Test normalization on missing variants
normalizer = VariantNormalizer('KCNH2')

# Read from missing_in_sqlite.csv
missing_file = 'comparison_results/missing_in_sqlite.csv'
total = 0
successful = 0

print('Testing normalization coverage...')

with open(missing_file, 'r') as f:
    reader = csv.DictReader(f)
    for i, row in enumerate(reader):
        if i >= 50:  # Test first 50
            break
        
        raw = row['excel_variant_raw'].strip()
        if not raw or 'Del(' in raw or 'Del(7)' in raw or 'ex' in raw.lower():
            continue  # Skip complex deletions for now
            
        normalized = normalizer.normalize_protein(raw)
        total += 1
        
        if normalized:
            successful += 1
            print(raw + ' -> ' + str(normalized))
        else:
            print(raw + ' -> FAILED')

print()
print('Coverage: ' + str(successful) + '/' + str(total) + ' = ' + str(int(100*successful/total)) + '%')