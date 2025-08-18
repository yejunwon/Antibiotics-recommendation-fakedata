import re
import json
import random
from docx import Document

# -------------------------------------
# 설정: 상품명 기준으로 명시적 매핑 테이블
# -------------------------------------
ingredient_to_product = {
    "Piperacillin/Tazobactam": "Tazoferan(R) 4.5g",
    "Cefazolin": "cefaZOLin 1g",
    "Azithromycin": "Azithromycin 250mg",
    "Ceftriaxone": "cefTRIAXone sod 2g",
    "Cefepime": "cefePIMe 1g",
    "Amoxicillin/Clavulanic Acid": "Amoxclan duo(R) 437.5mg/62.5mg",
    "Meropenem": "Meropenem 500mg"
}

# 상품명 리스트
target_products = list(ingredient_to_product.values())

# 감염균별 gram 분류
gram_mapping = {
    "Escherichia coli": "gram_negative",
    "Klebsiella pneumoniae": "gram_negative",
    "Klebsiella pneumoniae ssp pneumoniae": "gram_negative",
    "Enterobacter cloacae": "gram_negative",
    "Enterobacter cloacae ssp cloacae": "gram_negative",
    "Pseudomonas aeruginosa": "gram_negative",
    "Staphylococcus aureus": "gram_positive",
    "Enterococcus faecalis": "gram_positive"
}

severity_levels = ["mild", "moderate", "severe", "sepsis", "septic_shock"]

def sample_from_products(k=2):
    return random.sample(target_products, k=min(k, len(target_products)))

# -------------------------------------
# .docx 로드 및 유효 블록 파싱
# -------------------------------------
doc = Document("D:\dev_project\shl\Antibiotics\example.docx")
full_text = "\n".join([p.text for p in doc.paragraphs])
blocks = re.split(r"={10,}", full_text)

# 유효 감염균+감수성 결과만 저장
organism_sus_map = []

for block in blocks:
    org_match = re.search(r"Organism\s*:\s*(.+)", block)
    if not org_match:
        continue
    organism = org_match.group(1).strip()

    abx_lines = re.findall(r"^\s*(.+?)\s{2,}[\d<>=\.]+\s{2,}([SIR])\s*$", block, re.MULTILINE)
    mapped = {}

    for abx_name, sir in abx_lines:
        abx_name = abx_name.strip()
        if abx_name in ingredient_to_product:
            product_name = ingredient_to_product[abx_name]
            mapped[product_name] = sir

    if mapped:
        organism_sus_map.append((organism, mapped))

# -------------------------------------
# 각 감염균 블록당 N명 환자 생성
# -------------------------------------
N = 5  # block당 환자 수
patients = []
pid_counter = 1

for organism, susceptibility in organism_sus_map:
    for i in range(N):
        age = random.randint(18, 90)
        creatinine = round(random.uniform(0.6, 3.5), 1)
        bilirubin = round(random.uniform(0.4, 2.5), 1)
        ast = random.randint(15, 80)
        alt = random.randint(15, 80)
        neutropenia = random.choice([True, False])
        severity = random.choice(severity_levels)
        allergy = sample_from_products(k=random.randint(0, 2))
        drug_inter = sample_from_products(k=random.randint(0, 2))

        patient = {
            "patient_id": f"P{pid_counter:03d}",
            "age": age,
            "renal_function": {"creatinine": creatinine},
            "hepatic_function": {"bilirubin": bilirubin, "AST": ast, "ALT": alt},
            "neutropenia": neutropenia,
            "infection_severity": severity,
            "allergy": allergy,
            "drug_interactions": drug_inter,
            "infectious_agent": organism,
            "gram_status": gram_mapping.get(organism, "unknown"),
            "susceptibility": {
                organism: susceptibility
            }
        }
        patients.append(patient)
        pid_counter += 1

# -------------------------------------
# JSON 파일 저장
# -------------------------------------
with open("more_generated_patients.json", "w", encoding="utf-8") as f:
    json.dump(patients, f, indent=2, ensure_ascii=False)

print(f"✅ 총 {len(patients)}명 환자 JSON 생성 완료 → generated_patients.json")
