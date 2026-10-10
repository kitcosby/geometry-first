from sage.all import *
# Verify: s = d xor p for the 8 table states.
# Table: (code, d,p,m, s_exp)
# OII: (0,1,1) s=1 (+)
# IOI: (1,0,1) s=1 (+)
# OOO: (0,0,0) s=0 (-)
# IIO: (1,1,0) s=0 (-)
# OIO: (0,1,0) s=1 (+)
# IOO: (1,0,0) s=1 (+)
# OOI: (0,0,1) s=0 (-)
# III: (1,1,1) s=0 (-)

states = [
    ('OII',0,1,1,1),
    ('IOI',1,0,1,1),
    ('OOO',0,0,0,0),
    ('IIO',1,1,0,0),
    ('OIO',0,1,0,1),
    ('IOO',1,0,0,1),
    ('OOI',0,0,1,0),
    ('III',1,1,1,0),
]

print("=== s = d xor p verification ===")
all_ok = True
for code,d,p,m,s_exp in states:
    s_calc = d ^ p  # XOR
    ok = (s_calc == s_exp)
    if not ok: all_ok=False
    print(f"{code}: d={d},p={p} -> s_calc={s_calc}, s_exp={s_exp} [{'OK' if ok else 'FAIL'}]")

print()
print(f"s = d xor p holds for all 8: {all_ok}")
print()
# Even parity: d xor p xor s = 0?
print("=== Even parity (d xor p xor s = 0) ===")
for code,d,p,m,s_exp in states:
    par = d ^ p ^ s_exp
    print(f"{code}: parity={par} {'(even)' if par==0 else '(ODD!)'}")

print()
print("CONCLUSION:")
print("The 8 table states satisfy s=d xor p, hence d xor p xor s=0 (even parity).")
print("They lie in the 8-dim even-parity subspace of the 16-dim 4-bit tensor product.")
print("The 8 missing states (odd parity) are not in the table.")

with open('/root/parity.txt','w') as f:
    f.write(f"parity_verified={all_ok}\n")
print("WROTE /root/parity.txt")
