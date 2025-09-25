import numpy as np
import re, json
from typing import Any, Dict, List, Optional

# ====================== 데이터 구조 (원본) ======================
score2onehot = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 8: 7, 9: 8, 10: 9, 12: 10, 15: 11, 16: 12, 20: 13, 25: 14}
CATEGORY_MIN = 1
CATEGORY_MAX = 5
MAX_ONEHOT = 14

abx_nodes = [
    "Tazoferan(R) 4.5g", "cefaZOLin 1g", "Azithromycin 250mg", "cefTRIAXone sod 2g",
    "cefePIMe 1g", "Amoxclan duo(R) 437.5mg/62.5mg", "Meropenem 500mg"
]

abx_to_gram = {
    "Tazoferan(R) 4.5g": ["gram_positive", "gram_negative"],
    "cefaZOLin 1g": ["gram_positive"],
    "Azithromycin 250mg": ["gram_positive"],
    "cefTRIAXone sod 2g": ["gram_negative"],
    "cefePIMe 1g": ["gram_negative"],
    "Amoxclan duo(R) 437.5mg/62.5mg": ["gram_positive", "gram_negative"],
    "Meropenem 500mg": ["gram_positive", "gram_negative"]
}

abx_risk = {
    "Tazoferan(R) 4.5g": {"age": 3, "creatinine": 4, "hepatic": 2},
    "cefaZOLin 1g": {"age": 2, "creatinine": 3, "hepatic": 1},
    "Azithromycin 250mg": {"age": 1, "creatinine": 0, "hepatic": 3},
    "cefTRIAXone sod 2g": {"age": 2, "creatinine": 1, "hepatic": 3},
    "cefePIMe 1g": {"age": 3, "creatinine": 5, "hepatic": 1},
    "Amoxclan duo(R) 437.5mg/62.5mg": {"age": 1, "creatinine": 3, "hepatic": 4},
    "Meropenem 500mg": {"age": 2, "creatinine": 4, "hepatic": 2}
}

infectious_agents = [
    "Escherichia coli", "Staphylococcus aureus", "NTM(Nontuberculous Mycobacteria)",
    "Klebsiella pneumoniae", "Pseudomonas aeruginosa", "Enterococcus faecalis",
    "Enterococcus faecium"
]

agent_to_gram = {
    "Escherichia coli": "gram_negative",
    "Staphylococcus aureus": "gram_positive",
    "NTM(Nontuberculous Mycobacteria)": "gram_positive",
    "Klebsiella pneumoniae": "gram_negative",
    "Pseudomonas aeruginosa": "gram_negative",
    "Enterococcus faecalis": "gram_positive",
    "Enterococcus faecium": "gram_positive"
}

# === NEW: 병원체 Morphology & Metabolism 태그 ===
# morphology: {'coccus','coccus_clusters','coccus_chains','coccus_pairs','bacillus','coccobacillus','curved','filamentous'}
# metabolism: {'obligate_aerobe','facultative_anaerobe','obligate_anaerobe','microaerophile'}
agent_to_morphology = {
    "Escherichia coli": "bacillus",
    "Staphylococcus aureus": "coccus_clusters",
    "NTM(Nontuberculous Mycobacteria)": "bacillus",
    "Klebsiella pneumoniae": "bacillus",
    "Pseudomonas aeruginosa": "bacillus",
    "Enterococcus faecalis": "coccus_chains",
    "Enterococcus faecium": "coccus_chains"
}

agent_to_metabolism = {
    "Escherichia coli": "facultative_anaerobe",
    "Staphylococcus aureus": "facultative_anaerobe",
    "NTM(Nontuberculous Mycobacteria)": "obligate_aerobe",
    "Klebsiella pneumoniae": "facultative_anaerobe",
    "Pseudomonas aeruginosa": "obligate_aerobe",
    "Enterococcus faecalis": "facultative_anaerobe",
    "Enterococcus faecium": "facultative_anaerobe"
}

# 새로 추가: 항생제 그룹화 (Cephalosporin 계열 Group A, 나머지 Group B)
abx_to_group = {
    "cefaZOLin 1g": "A",  # Cephalosporin
    "cefTRIAXone sod 2g": "A",  # Cephalosporin
    "cefePIMe 1g": "A",  # Cephalosporin
    "Tazoferan(R) 4.5g": "B",
    "Azithromycin 250mg": "B",
    "Amoxclan duo(R) 437.5mg/62.5mg": "B",
    "Meropenem 500mg": "B"
}

# 새로 추가: 그룹별 CrCl 값에 따른 creat_score 매핑 함수 (임의 값으로 예시)
def get_creat_score(crcl, group):
    """CrCl 값에 따라 그룹별 creat_score (1~5) 반환. 높은 score = 높은 위험."""
    if group == "A":  # Cephalosporin 계열: 더 엄격한 기준
        if crcl >= 90:
            return 1
        elif crcl >= 60:
            return 2
        elif crcl >= 30:  # (버그 수정: 3 -> 30)
            return 3
        elif crcl >= 15:
            return 4
        else:
            return 5
    elif group == "B":  # 나머지: 더 완화된 기준
        if crcl >= 100:
            return 1
        elif crcl >= 70:
            return 2
        elif crcl >= 40:
            return 3
        elif crcl >= 20:
            return 4
        else:
            return 5
    else:
        raise ValueError(f"Unknown group: {group}")

# Min-Max 정규화 함수 (기존 유지)
def normalize_score(score, min_val=0, max_val=MAX_ONEHOT):
    if score == 0:
        return 0.0
    if max_val <= min_val:
        print(f"Warning: max_val({max_val}) <= min_val({min_val}). Returning 0.")
        return 0
    return (score - min_val) / (max_val - min_val)

# 가중합 순위화 함수 (변경: creat_score를 abx 그룹에 따라 동적으로 계산)
def rank_by_weighted_sum(candidates, age_score, hepatic_score, crcl, weights=[0.2857, 0.4286, 0.2857]):
    if not candidates:
        return [], [], 0.0
    
    # 가중치 검증
    if abs(sum(weights) - 1.0) > 1e-6:
        print(f"Warning: Sum of weights({sum(weights)}) != 1. Normalizing.")
        weights = [w / sum(weights) for w in weights]
    
    # 각 abx별로 creat_score 동적으로 계산
    data = np.array([
        [
            normalize_score(score2onehot.get(age_score * abx_risk[abx]['age'], 0)),
            normalize_score(score2onehot.get(get_creat_score(crcl, abx_to_group[abx]) * abx_risk[abx]['creatinine'], 0)),
            normalize_score(score2onehot.get(hepatic_score * abx_risk[abx]['hepatic'], 0))
        ]
        for abx in candidates
    ])
    
    # 가중합 계산
    weighted_scores = np.dot(data, np.array(weights))
    total_score = np.sum(weighted_scores)
    sorted_idx = np.argsort(weighted_scores)  # 낮은 점수 순 (부작용 적음)
    
    # 로그
    log = [
        f"{candidates[i]} (Score={weighted_scores[i]:.3f}, weighted_age={data[i][0] * weights[0]:.3f}, weighted_creat={data[i][1] * weights[1]:.3f}, weighted_hepatic={data[i][2] * weights[2]:.3f})"
        for i in sorted_idx
    ]
    
    return [candidates[i] for i in sorted_idx], log, total_score

# 환자 상태 요약 함수 (변경: extreme_creatinine 제거, 그룹별 다르기 때문)
def get_patient_states(age_score, hepatic_score):
    states = []
    if age_score >= 4:
        states.append("extreme_age")
    if hepatic_score >= 4:
        states.append("extreme_hepatic")
    return states

# 상태 점수 함수 (변경: creat_score 제거)
def get_status_score(age_category, hepatic_category):
    return (
        validate_category(age_category, "나이 카테고리"),
        validate_category(hepatic_category, "빌리루빈 카테고리")
    )

# 입력 검증 함수 (변경: CrCl 숫자 검증 추가)
def validate_category(value, field_name):
    try:
        val = int(value)
        if CATEGORY_MIN <= val <= CATEGORY_MAX:
            return val
        raise ValueError(f"{field_name}는 {CATEGORY_MIN}~{CATEGORY_MAX} 사이여야 합니다.")
    except ValueError:
        raise ValueError(f"{field_name}는 숫자여야 하며 {CATEGORY_MIN}~{CATEGORY_MAX} 사이여야 합니다.")

def validate_crcl(value):
    try:
        val = float(value)
        if val > 0:
            return val
        raise ValueError("CrCl 값은 0보다 커야 합니다.")
    except ValueError:
        raise ValueError("CrCl 값은 숫자여야 합니다.")

def validate_antibiotics(input_list, valid_list, field_name):
    valid_items = [item.strip() for item in input_list if item.strip() in valid_list]
    if len(valid_items) < len([item for item in input_list if item.strip()]):
        print(f"경고: 일부 {field_name}가 유효하지 않습니다. 유효한 항목만 포함됩니다.")
    return valid_items

# ====================== NEW: Taxonomy RegEx 탐지기 ======================
# 형태학 키워드
_RE_GPC = re.compile(r'\b(?:gpc|gram\s*positive\s*cocci|g\+\s*cocci)\b', re.I)
_RE_GPR = re.compile(r'\b(?:gpr|gram\s*positive\s*rod[s]?|g\+\s*rod[s]?|gram\s*positive\s*bacill[i]?)\b', re.I)
_RE_GNR = re.compile(r'\b(?:gnr|gram\s*negative\s*rod[s]?|g-\s*rod[s]?|gram\s*negative\s*bacill[i]?)\b', re.I)
_RE_GNC = re.compile(r'\b(?:gnc|gram\s*negative\s*cocci)\b', re.I)

_RE_COCCI = re.compile(r'\b(cocci|coccus)\b', re.I)
_RE_ROD   = re.compile(r'\b(rod[s]?|bacill[i]?)\b', re.I)
_RE_COCCOBAC = re.compile(r'\b(coccobacill[i]?)\b', re.I)
_RE_CURVED = re.compile(r'\b(curved|comma[-\s]*shaped|vibrio|campylobacter)\b', re.I)
_RE_FILAMENT = re.compile(r'\b(filamentous|branching|actinomyces|nocardia)\b', re.I)

_RE_CLUSTER = re.compile(r'\b(cluster[s]?|in\s*clusters)\b', re.I)
_RE_CHAIN   = re.compile(r'\b(chain[s]?|in\s*chains)\b', re.I)
_RE_PAIRS   = re.compile(r'\b(pair[s]?|in\s*pairs|diplococci)\b', re.I)

# 대사 키워드
_RE_OBL_AER = re.compile(r'\bobligate\s*aerobe[s]?\b', re.I)
_RE_FAC_ANA = re.compile(r'\bfacultative\s*anaerobe[s]?\b', re.I)
_RE_OBL_ANA = re.compile(r'\bobligate\s*anaerobe[s]?\b', re.I)
_RE_MICRO   = re.compile(r'\bmicroaerophil\w*\b', re.I)
_RE_AEROBIC = re.compile(r'\baerobic\b', re.I)
_RE_ANAEROB = re.compile(r'\banaerob(ic|e)s?\b', re.I)

def _detect_taxonomy_from_text(text: str) -> Dict[str, Optional[str]]:
    """Free-text에서 Gram/Morphology/Metabolism을 정규식으로 추정."""
    gram: Optional[str] = None
    morph: Optional[str] = None
    metab: Optional[str] = None

    # Gram + Morph combo shorthands
    if _RE_GPC.search(text):
        gram = "gram_positive"
        morph = "coccus"
    if _RE_GPR.search(text):
        gram = "gram_positive"
        morph = "bacillus"
    if _RE_GNR.search(text):
        gram = "gram_negative"
        morph = "bacillus"
    if _RE_GNC.search(text):
        gram = "gram_negative"
        morph = "coccus"

    # Fallback morphology words
    if morph is None:
        if _RE_COCCOBAC.search(text):
            morph = "coccobacillus"
        elif _RE_ROD.search(text):
            morph = "bacillus"
        elif _RE_COCCI.search(text):
            morph = "coccus"
        elif _RE_CURVED.search(text):
            morph = "curved"
        elif _RE_FILAMENT.search(text):
            morph = "filamentous"

    # Arrangement refinement for cocci
    if morph and "coccus" in morph:
        if _RE_CLUSTER.search(text):
            morph = "coccus_clusters"
        elif _RE_CHAIN.search(text):
            morph = "coccus_chains"
        elif _RE_PAIRS.search(text):
            morph = "coccus_pairs"

    # Metabolism
    if _RE_OBL_AER.search(text):
        metab = "obligate_aerobe"
    elif _RE_FAC_ANA.search(text):
        metab = "facultative_anaerobe"
    elif _RE_OBL_ANA.search(text):
        metab = "obligate_anaerobe"
    elif _RE_MICRO.search(text):
        metab = "microaerophile"
    else:
        # generic words
        if _RE_AEROBIC.search(text):
            metab = "obligate_aerobe"
        elif _RE_ANAEROB.search(text):
            metab = "obligate_anaerobe"

    return {"gram_status": gram, "morphology": morph, "metabolism": metab}

# 항생제 추천 함수 (감염균 알려진 경우, 변경: creat_score 동적, extreme_creatinine abx별 판단)
def recommend_antibiotics(patient):
    gram = patient.get('gram_status')
    agent = patient.get('infectious_agent')
    allergy = set(patient.get('allergy', []))
    drug_inter = set(patient.get('drug_interactions', []))
    suscept = patient.get('susceptibility', {}).get(agent, {})
    crcl = patient['renal_function']['creatinine_value']
    log = ["━━━━━━━━━━━━━━━━━━━━━━━", f"🔍 환자 상태 요약: {agent if agent else 'Unknown agent'} ({gram if gram else 'unknown'})"]

    # NEW: morphology/metabolism 표시
    morph = patient.get('morphology', 'unknown')
    metab = patient.get('metabolism', 'unknown')
    log.append(f"  - Morphology/Metabolism: {morph} / {metab}")

    # 환자 상태 (creat_score 없음, 동적)
    age_score, hepatic_score = get_status_score(
        patient['age_category'], 
        patient['hepatic_function']['bilirubin_category']
    )
    patient_states = get_patient_states(age_score, hepatic_score)
    state_str = ", ".join(patient_states) if patient_states else "없음"
    log.append(f"  - 극도의 상태 (공통): {state_str}")
    log.append(f"  - CrCl 값: {crcl}")
    log.append(f"  - 알러지: {', '.join(allergy) if allergy else '없음'}")
    log.append(f"  - 약물 상호작용: {', '.join(drug_inter) if drug_inter else '없음'}")
    log.append(f"  - 감수성: {', '.join([f'{k}: {v}' for k, v in suscept.items()]) if suscept else '없음'}")
    log.append(f"  - 가중치: 나이={0.2857:.4f}, 신장={0.4286:.4f}, 간={0.2857:.4f}")

    # 1. Gram 상태 필터링 (기존)
    candidates = [abx for abx in abx_nodes if gram in abx_to_gram[abx]] if gram else list(abx_nodes)
    eliminated_gram = [f"  · {abx}: {gram}에 효과 없음" for abx in abx_nodes if gram and abx not in candidates]
    log.append("🔹 1. Gram 상태 필터링")
    log.append(f"  - 후보: {', '.join(candidates) if candidates else '없음'}")
    log.append("  ⮩ [제외]: " + (", ".join(eliminated_gram) if eliminated_gram else "없음"))

    # 2. 독성 필터링 (변경: abx별 creat_score로 판단)
    filtered = []
    eliminated_state = []
    for abx in candidates:
        creat_score = get_creat_score(crcl, abx_to_group[abx])
        toxic_reasons = []
        if "extreme_age" in patient_states and abx_risk[abx]['age'] >= 4:
            toxic_reasons.append("extreme_age")
        if creat_score >= 4 and abx_risk[abx]['creatinine'] >= 4:
            toxic_reasons.append("extreme_creatinine")
        if "extreme_hepatic" in patient_states and abx_risk[abx]['hepatic'] >= 4:
            toxic_reasons.append("extreme_hepatic")
        if toxic_reasons:
            eliminated_state.append(f"  · {abx}: 독성 위험 ({', '.join(toxic_reasons)}), creat_score={creat_score}")
        else:
            filtered.append(abx)
    log.append("🔹 2. 독성 필터링")
    log.append(f"  - 후보: {', '.join(filtered) if filtered else '없음'}")
    log.append("  ⮩ [제외]: " + (", ".join(eliminated_state) if eliminated_state else "없음"))

    # 3. 알러지/약물 상호작용 필터링 (기존)
    filtered2 = []
    eliminated_allergy = []
    for abx in filtered:
        reasons = []
        if abx in allergy:
            reasons.append("알러지")
        if abx in drug_inter:
            reasons.append("약물 상호작용")
        if reasons:
            eliminated_allergy.append(f"  · {abx}: {', '.join(reasons)}")
        else:
            filtered2.append(abx)
    log.append("🔹 3. 알러지/상호작용 필터링")
    log.append(f"  - 최종 후보: {', '.join(filtered2) if filtered2 else '없음'}")
    log.append("  ⮩ [제외]: " + (", ".join(eliminated_allergy) if eliminated_allergy else "없음"))

    # 4. 감수성 필터링 (기존)
    final_candidates = []
    eliminated_suscept = []
    for abx in filtered2:
        if suscept.get(abx, None) == "S":
            final_candidates.append(abx)
        else:
            eliminated_suscept.append(f"  · {abx}: 감수성 미달 ({suscept.get(abx, '없음')})")
    if not final_candidates:
        for abx in filtered2:
            if suscept.get(abx, None) == "I":
                final_candidates.append(abx)
                log.append(f"  · {abx}: 감수성 I (S 없으므로 포함)")
            else:
                eliminated_suscept.append(f"  · {abx}: 감수성 미달 ({suscept.get(abx, '없음')})")
    log.append("🔹 4. 감수성 필터링")
    log.append(f"  - 최종 후보: {', '.join(final_candidates) if final_candidates else '없음'}")
    log.append("  ⮩ [제외]: " + (", ".join(eliminated_suscept) if eliminated_suscept else "없음"))

    # 5. 가중합 부작용 순위 (변경: crcl 전달)
    if final_candidates:
        log.append("🔹 5. 가중합 부작용 순위")
        ranked_candidates, weighted_log, total_score = rank_by_weighted_sum(
            final_candidates, age_score, hepatic_score, crcl
        )
        log.append("  - 순위:")
        log.extend([f"    · {res}" for res in weighted_log])
    else:
        log.append("⚠ 추천 항생제 없음")
        ranked_candidates = []
        total_score = 0.0

    log.append("━━━━━━━━━━━━━━━━━━━━━━━")
    return ranked_candidates, log

# Empiric 항생제 추천 함수 (변경 유사)
def recommend_empiric_antibiotics(patient):
    log = ["━━━━━━━━━━━━━━━━━━━━━━━", "🔍 Empiric 항생제 추천 (감염균 미확인)"]
    age_score, hepatic_score = get_status_score(
        patient['age_category'], 
        patient['hepatic_function']['bilirubin_category']
    )
    crcl = patient['renal_function']['creatinine_value']
    patient_states = get_patient_states(age_score, hepatic_score)
    allergy = set(patient.get('allergy', []))
    drug_inter = set(patient.get('drug_interactions', []))

    # NEW: morphology/metabolism 표시
    morph = patient.get('morphology', 'unknown')
    metab = patient.get('metabolism', 'unknown')
    gram = patient.get('gram_status', 'unknown')
    log.append(f"  - Gram/Morphology/Metabolism: {gram} / {morph} / {metab}")

    state_str = ", ".join(patient_states) if patient_states else "없음"
    log.append(f"  - 극도의 상태 (공통): {state_str}")
    log.append(f"  - CrCl 값: {crcl}")
    log.append(f"  - 알러지: {', '.join(allergy) if allergy else '없음'}")
    log.append(f"  - 약물 상호작용: {', '.join(drug_inter) if drug_inter else '없음'}")
    log.append(f"  - 가중치: 나이={0.2857:.4f}, 신장={0.4286:.4f}, 간={0.2857:.4f}")

    # 1. 초기 후보 (기존)
    candidates = [abx for abx in abx_nodes if len(abx_to_gram[abx]) == 2]
    candidates += [abx for abx in abx_nodes if abx not in candidates]
    log.append("🔹 1. 초기 후보")
    log.append(f"  - 후보: {', '.join(candidates)}")

    # 2. 독성 필터링 (변경: abx별 creat_score)
    filtered = []
    eliminated_state = []
    for abx in candidates:
        creat_score = get_creat_score(crcl, abx_to_group[abx])
        toxic_reasons = []
        if "extreme_age" in patient_states and abx_risk[abx]['age'] >= 4:
            toxic_reasons.append("extreme_age")
        if creat_score >= 4 and abx_risk[abx]['creatinine'] >= 4:
            toxic_reasons.append("extreme_creatinine")
        if "extreme_hepatic" in patient_states and abx_risk[abx]['hepatic'] >= 4:
            toxic_reasons.append("extreme_hepatic")
        if toxic_reasons:
            eliminated_state.append(f"  · {abx}: 독성 위험 ({', '.join(toxic_reasons)}), creat_score={creat_score}")
        else:
            filtered.append(abx)
    log.append("🔹 2. 독성 필터링")
    log.append(f"  - 후보: {', '.join(filtered) if filtered else '없음'}")
    log.append("  ⮩ [제외]: " + (", ".join(eliminated_state) if eliminated_state else "없음"))

    # 3. 알러지/약물 상호작용 필터링 (기존)
    final_candidates = []
    eliminated_allergy = []
    for abx in filtered:
        reasons = []
        if abx in allergy:
            reasons.append("알러지")
        if abx in drug_inter:
            reasons.append("약물 상호작용")
        if reasons:
            eliminated_allergy.append(f"  · {abx}: {', '.join(reasons)}")
        else:
            final_candidates.append(abx)
    log.append("🔹 3. 알러지/상호작용 필터링")
    log.append(f"  - 최종 후보: {', '.join(final_candidates) if final_candidates else '없음'}")
    log.append("  ⮩ [제외]: " + (", ".join(eliminated_allergy) if eliminated_allergy else "없음"))

    # 4. 가중합 부작용 순위 (변경: crcl 전달)
    if final_candidates:
        log.append("🔹 4. 가중합 부작용 순위")
        ranked_candidates, weighted_log, total_score = rank_by_weighted_sum(
            final_candidates, age_score, hepatic_score, crcl
        )
        log.append("  - 순위:")
        log.extend([f"    · {res}" for res in weighted_log])
    else:
        log.append("⚠ 추천 항생제 없음")
        ranked_candidates = []
        total_score = 0.0

    log.append("━━━━━━━━━━━━━━━━━━━━━━━")
    return ranked_candidates, log

# 사용자 입력 함수 (원본: 이제 Paste-&-Go에서 대체 사용, 필요시 유지)
def get_user_input():
    print("=== 환자 데이터 입력 ===")
    pathogen_known = input("감염균을 아나요? (y/n): ").strip().lower() == 'y'

    print("나이 카테고리 (1~5): 1: 10세 이하, 2: 11~30세, 3: 31~50세, 4: 51~60세, 5: 61세 이상")
    age_category = validate_category(input("선택: "), "나이 카테고리")

    print("\nCrCl 값 (mL/min, 숫자 입력, 예: 85.5):")
    creatinine_value = validate_crcl(input("입력: "))

    print("\n빌리루빈 카테고리 (1~5): 1: 0.9 이하, 2: 0.9~1.2, 3: 1.2~1.8, 4: 1.8~2.3, 5: 2.3 초과")
    bilirubin_category = validate_category(input("선택: "), "빌리루빈 카테고리")

    print("\n알러지 항생제 (쉼표로 구분, 없으면 엔터, 유효 항생제: " + ", ".join(abx_nodes) + "):")
    allergy_input = input().strip()
    allergy = validate_antibiotics([a.strip() for a in allergy_input.split(",") if a.strip()], abx_nodes, "알러지")

    print("\n약물 상호작용 항생제 (쉼표로 구분, 없으면 엔터, 유효 항생제: " + ", ".join(abx_nodes) + "):")
    drug_inter_input = input().strip()
    drug_inter = validate_antibiotics([d.strip() for d in drug_inter_input.split(",") if d.strip()], abx_nodes, "약물 상호작용")

    patient = {
        "age_category": age_category,
        "renal_function": {"creatinine_value": creatinine_value},
        "hepatic_function": {"bilirubin_category": bilirubin_category},
        "allergy": allergy,
        "drug_interactions": drug_inter,
        "pathogen_known": pathogen_known
    }

    if pathogen_known:
        print("\n감염균 선택:")
        for i, agent in enumerate(infectious_agents, 1):
            print(f"{i}. {agent}")
        while True:
            try:
                agent_idx = int(input(f"선택 (1~{len(infectious_agents)}): "))
                if 1 <= agent_idx <= len(infectious_agents):
                    break
                print(f"1~{len(infectious_agents)} 사이의 숫자를 입력하세요.")
            except ValueError:
                print("숫자를 입력하세요.")
        infectious_agent = infectious_agents[agent_idx - 1]
        patient["infectious_agent"] = infectious_agent
        patient["gram_status"] = agent_to_gram[infectious_agent]
        # NEW: morphology/metabolism 태깅
        patient["morphology"] = agent_to_morphology.get(infectious_agent)
        patient["metabolism"] = agent_to_metabolism.get(infectious_agent)

        susceptibility = {}
        print("\n항생제별 감수성(S, I, R, 기본값 S):")
        for abx in abx_nodes:
            sir = input(f"{abx}: ").strip().upper()
            susceptibility[abx] = "S" if sir == "" else sir
            if sir not in ["", "S", "I", "R"]:
                print(f"경고: {abx}의 감수성 {sir}은 유효하지 않음. 기본값 S로 설정.")
        patient["susceptibility"] = {infectious_agent: susceptibility}

    return patient

# ====================== Paste-&-Go 입력부 (붙여넣기 파서) ======================
# RegEx 패턴 (한/영 혼용, 느슨매칭)
_RE_AGE_NUM = re.compile(r'\b(?:age|나이)\s*[:=]?\s*(\d{1,3})\b', re.I)
_RE_AGE_YO  = re.compile(r'\b(\d{1,3})\s*(?:yo|y/o|yr|yrs|years?\s*old)\b', re.I)
_RE_CRCL    = re.compile(r'\b(?:CrCl|creatinine\s*clearance|크레아티닌\s*청소율)\s*[:=]?\s*([0-9]+(?:\.[0-9]+)?)\s*(?:mL/?min|ml/?min)?', re.I)
_RE_BILI    = re.compile(r'\b(?:total\s*bilirubin|bilirubin|TBili|T\.?bil|총\s*빌리루빈)\s*[:=]?\s*([0-9]+(?:\.[0-9]+)?)\s*(?:mg/dL|mg\s*/\s*dL)?', re.I)
_RE_SUSC_LINE = re.compile(r'^\s*([A-Za-z][^:]+?)\s*:\s*([SIR])\s*$', re.I | re.M)
_RE_ALLERGY = re.compile(r'(?:^|\n)\s*(?:allerg(?:y|ies)|알러지)\s*[:=]\s*(.+)', re.I)
_RE_DRUGINT = re.compile(r'(?:^|\n)\s*(?:drug\s*interactions?|상호작용)\s*[:=]\s*(.+)', re.I)

def _scan_abx_names(fragment: str) -> List[str]:
    hits, low = [], fragment.lower()
    for abx in abx_nodes:
        if abx.lower() in low:
            hits.append(abx)
    return list(dict.fromkeys(hits))

def _detect_agent(text: str) -> Optional[str]:
    low = text.lower()
    for ag in infectious_agents:
        if ag.lower() in low:
            return ag
    return None

def _age_to_category(age_years: int) -> int:
    if age_years <= 10: return 1
    if age_years <= 30: return 2
    if age_years <= 50: return 3
    if age_years <= 60: return 4
    return 5

def _bili_to_category(tbili: float) -> int:
    if tbili <= 0.9: return 1
    if tbili <= 1.2: return 2
    if tbili <= 1.8: return 3
    if tbili <= 2.3: return 4
    return 5

# NEW: 텍스트에서 Gram/Morph/Metabolism 추출
def _detect_taxonomy(text: str) -> Dict[str, Optional[str]]:
    return _detect_taxonomy_from_text(text)

def _extract_from_text(emr_text: str) -> Dict[str, Any]:
    text = emr_text.strip()

    # age
    age_val = None
    m = _RE_AGE_NUM.search(text) or _RE_AGE_YO.search(text)
    if m:
        gs = [g for g in m.groups() if g]
        if gs:
            age_val = int(gs[0])

    # CrCl
    crcl_val = None
    m = _RE_CRCL.search(text)
    if m:
        crcl_val = float(m.group(1))

    # Total bilirubin
    tbili_val = None
    m = _RE_BILI.search(text)
    if m:
        tbili_val = float(m.group(1))

    # 알러지/상호작용
    allergy: List[str] = []
    m = _RE_ALLERGY.search(text)
    if m:
        allergy = _scan_abx_names(m.group(1))

    drug_inter: List[str] = []
    m = _RE_DRUGINT.search(text)
    if m:
        drug_inter = _scan_abx_names(m.group(1))

    # 감염균/Gram/감수성
    infectious_agent: Optional[str] = _detect_agent(text)

    # NEW: 우선 텍스트에서 Gram/Morph/Metabolism 추정
    taxo = _detect_taxonomy(text)
    gram_status: Optional[str] = taxo.get("gram_status")
    morphology: Optional[str] = taxo.get("morphology")
    metabolism: Optional[str] = taxo.get("metabolism")

    # 병원체로 override/보완
    if infectious_agent:
        gram_status = agent_to_gram[infectious_agent]
        morphology = agent_to_morphology.get(infectious_agent, morphology)
        metabolism = agent_to_metabolism.get(infectious_agent, metabolism)

    susceptibility = {}
    if infectious_agent:
        susc_map = {}
        for m in _RE_SUSC_LINE.finditer(text):
            abx_name = m.group(1).strip()
            sir = m.group(2).upper()
            for abx in abx_nodes:
                if abx_name.lower() in abx.lower() or abx.lower() in abx_name.lower():
                    if sir in {"S","I","R"}:
                        susc_map[abx] = sir
        if susc_map:
            susceptibility = {infectious_agent: susc_map}

    # 환자 딕트
    patient: Dict[str, Any] = {
        "allergy": allergy or [],
        "drug_interactions": drug_inter or [],
    }

    if age_val is not None:
        patient["age_category"] = _age_to_category(age_val)
    if tbili_val is not None:
        patient["hepatic_function"] = {"bilirubin_category": _bili_to_category(tbili_val)}
    if crcl_val is not None:
        patient["renal_function"] = {"creatinine_value": crcl_val}

    # NEW: taxonomy 태깅 결과 반영
    if gram_status:
        patient["gram_status"] = gram_status
    if morphology:
        patient["morphology"] = morphology
    if metabolism:
        patient["metabolism"] = metabolism

    if infectious_agent:
        patient["pathogen_known"] = True
        patient["infectious_agent"] = infectious_agent
        # susceptibility
        if susceptibility:
            patient["susceptibility"] = susceptibility
    else:
        patient["pathogen_known"] = False

    return patient

def _extract_from_json_str(emr_json_str: str) -> Optional[Dict[str, Any]]:
    try:
        data = json.loads(emr_json_str.strip())
    except Exception:
        return None
    if not isinstance(data, (dict, list)):
        return None
    if isinstance(data, list):
        data = next((x for x in data if isinstance(x, dict)), None)
        if data is None:
            return None

    patient: Dict[str, Any] = {
        "allergy": data.get("allergy", []) or [],
        "drug_interactions": data.get("drug_interactions", []) or [],
    }

    # age
    age_val = data.get("age") or data.get("Age")
    if isinstance(age_val, (int, float)):
        patient["age_category"] = _age_to_category(int(age_val))

    # CrCl
    crcl = None
    renal = data.get("renal") or data.get("kidney") or {}
    if isinstance(renal, dict):
        crcl = renal.get("CrCl") or renal.get("crcl") or renal.get("creatinine_clearance")
    if crcl is not None:
        patient["renal_function"] = {"creatinine_value": float(crcl)}

    # bilirubin
    tbili = None
    biomarkers = data.get("biomarkers") or {}
    if isinstance(biomarkers, dict):
        tbili = biomarkers.get("bilirubin") or biomarkers.get("TBili") or biomarkers.get("total_bilirubin")
    if tbili is not None:
        patient["hepatic_function"] = {"bilirubin_category": _bili_to_category(float(tbili))}

    # pathogen & susceptibility
    pathogen = data.get("pathogen")
    # NEW: taxonomy from free text field (optional)
    free_text = data.get("note") or data.get("text") or ""
    taxo = _detect_taxonomy(free_text) if isinstance(free_text, str) else {}

    if pathogen in agent_to_gram:
        patient["pathogen_known"] = True
        patient["infectious_agent"] = pathogen
        patient["gram_status"] = agent_to_gram[pathogen]
        patient["morphology"] = agent_to_morphology.get(pathogen, taxo.get("morphology"))
        patient["metabolism"] = agent_to_metabolism.get(pathogen, taxo.get("metabolism"))

        susc_in = data.get("susceptibility") or {}
        norm = {}
        if isinstance(susc_in, dict):
            for k, v in susc_in.items():
                v2 = str(v).upper()
                if v2 not in {"S","I","R"}:
                    continue
                for abx in abx_nodes:
                    if k.lower() in abx.lower() or abx.lower() in k.lower():
                        norm[abx] = v2
        if norm:
            patient["susceptibility"] = {pathogen: norm}
    else:
        patient["pathogen_known"] = False
        # 병원체 미확정이면 RegEx taxonomy 반영 시도
        if taxo.get("gram_status"):
            patient["gram_status"] = taxo.get("gram_status")
        if taxo.get("morphology"):
            patient["morphology"] = taxo.get("morphology")
        if taxo.get("metabolism"):
            patient["metabolism"] = taxo.get("metabolism")

    return patient

def build_patient_from_paste(raw_input: str) -> Dict[str, Any]:
    """자연어 또는 JSON을 그대로 붙여넣으면 파싱해서 patient dict 구성.
       필수값 누락 시 해당 항목만 추가 질문."""
    raw = raw_input.strip()
    patient = _extract_from_json_str(raw) or _extract_from_text(raw)

    # 기본 보정
    patient.setdefault("allergy", [])
    patient.setdefault("drug_interactions", [])
    if "renal_function" in patient:
        patient["renal_function"].setdefault("creatinine_value", None)

    # 누락 체크
    missing = []
    if "age_category" not in patient:
        missing.append("age_category")
    if not patient.get("renal_function", {}).get("creatinine_value"):
        missing.append("crcl")
    if "hepatic_function" not in patient or "bilirubin_category" not in patient["hepatic_function"]:
        missing.append("bilirubin_category")

    # 필요한 것만 추가 질의
    if missing:
        print("⚠ 일부 정보가 부족합니다. 아래 항목만 추가 입력을 받습니다.")
        if "age_category" in missing:
            print("나이 카테고리 (1~5): 1: ≤10, 2: 11–30, 3: 31–50, 4: 51–60, 5: ≥61")
            while True:
                try:
                    val = int(input("선택: ").strip())
                    validate_category(val, "나이 카테고리")
                    patient["age_category"] = val
                    break
                except Exception as e:
                    print(e)
        if "crcl" in missing:
            while True:
                try:
                    val = float(input("CrCl (mL/min): ").strip())
                    validate_crcl(val)
                    patient["renal_function"] = {"creatinine_value": val}
                    break
                except Exception as e:
                    print(e)
        if "bilirubin_category" in missing:
            print("빌리루빈 카테고리 (1~5): 1: ≤0.9, 2: 0.9–1.2, 3: 1.2–1.8, 4: 1.8–2.3, 5: >2.3")
            while True:
                try:
                    val = int(input("선택: ").strip())
                    validate_category(val, "빌리루빈 카테고리")
                    patient["hepatic_function"] = {"bilirubin_category": val}
                    break
                except Exception as e:
                    print(e)

    # 감염균은 등록 리스트에 있을 때만 인정. 없으면 empiric
    if not patient.get("pathogen_known"):
        patient["pathogen_known"] = False

    return patient

# ====================== 실행부 (Paste-&-Go 메인) ======================
MOCK_EMR_TEXT_COMPLETE = """
Age: 67
CrCl: 48 mL/min
Total bilirubin: 1.4 mg/dL
Allergy: cefaZOLin 1g
Drug interactions: Meropenem 500mg
Findings: GNR with foul-smelling abscess (likely anaerobe)
Pathogen: Escherichia coli
cefTRIAXone: S
cefePIMe 1g: R
Amoxclan duo: I
""".strip()

MOCK_EMR_JSON_COMPLETE = json.dumps({
    "age": 62,
    "renal": {"CrCl": 55},
    "biomarkers": {"bilirubin": 1.1},
    "allergy": ["Tazoferan(R) 4.5g"],
    "drug_interactions": [],
    "note": "GPC in clusters suspected; aerobic",
    "pathogen": "Klebsiella pneumoniae",
    "susceptibility": {"cefTRIAXone": "S", "cefePIMe": "I"}
}, ensure_ascii=False, indent=2)

def _run_once(raw: str):
    patient = build_patient_from_paste(raw)
    print("\n[PARSED PATIENT]\n", patient)
    print("\n=== 추천 결과 ===")
    if patient.get('pathogen_known'):
        result, log = recommend_antibiotics(patient)
    else:
        result, log = recommend_empiric_antibiotics(patient)
    print("\n".join(log))
    print("\n추천 항생제:", ", ".join(result) if result else "없음")

def main():
    print("=== EMR를 그대로 붙여넣으세요 (자연어 또는 JSON).")
    print("여러 명을 한 번에 처리하려면 환자 간 구분선으로 '---' 한 줄을 넣으세요.")
    print("붙여넣기 종료는 빈 줄 2번(연속)입니다. ===\n")

    while True:
        # 입력 수집
        print("\nPaste and enter twice")
        lines = []
        blank_streak = 0
        try:
            while True:
                line = input()
                if line.strip() == "":
                    blank_streak += 1
                    if blank_streak >= 2:
                        break
                else:
                    blank_streak = 0
                lines.append(line)
        except EOFError:
            pass

        raw_all = "\n".join(lines).strip()

        # 아무 것도 안 넣었으면 데모 실행
        if not raw_all:
            print("\n[DEMO MODE] → 텍스트/JSON 예제로 2건 실행")
            for demo in (MOCK_EMR_TEXT_COMPLETE, MOCK_EMR_JSON_COMPLETE):
                print("\n[DEMO INPUT]\n" + demo)
                _run_once(demo)
        else:
            # 배치 지원: '---' (단독 라인)로 나눈다
            chunks = []
            buf = []
            for ln in raw_all.splitlines():
                if ln.strip() == "---":
                    if buf:
                        chunks.append("\n".join(buf).strip())
                        buf = []
                else:
                    buf.append(ln)
            if buf:
                chunks.append("\n".join(buf).strip())

            # 각 환자 처리
            for idx, chunk in enumerate(chunks, 1):
                print(f"\n===== 환자 {idx} =====")
                _run_once(chunk)

        # 계속 여부
        ans = input("\n계속 입력하시겠습니까? (y/n): ").strip().lower()
        if ans != "y":
            print("종료합니다.")
            break


if __name__ == "__main__":
    main()
