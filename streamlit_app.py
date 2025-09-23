import numpy as np

# One-hot encoding 매핑 (0~25 -> 0~14)
score2onehot = {
    0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 8: 7, 9: 8,
    10: 9, 12: 10, 15: 11, 16: 12, 20: 13, 25: 14
}

# 상수 정의
CATEGORY_MIN = 1
CATEGORY_MAX = 5
MAX_ONEHOT = 14

def get_onehot(score):
    return score2onehot.get(score, 0)

# 항생제 목록
abx_nodes = [
    "Tazoferan(R) 4.5g",
    "cefaZOLin 1g",
    "Azithromycin 250mg",
    "cefTRIAXone sod 2g",
    "cefePIMe 1g",
    "Amoxclan duo(R) 437.5mg/62.5mg",
    "Meropenem 500mg"
]

# 그람 양성/음성 적용 범위
abx_to_gram = {
    "Tazoferan(R) 4.5g": ["gram_positive", "gram_negative"],
    "cefaZOLin 1g": ["gram_positive"],
    "Azithromycin 250mg": ["gram_positive"],
    "cefTRIAXone sod 2g": ["gram_negative"],
    "cefePIMe 1g": ["gram_negative"],
    "Amoxclan duo(R) 437.5mg/62.5mg": ["gram_positive", "gram_negative"],
    "Meropenem 500mg": ["gram_positive", "gram_negative"]
}

# 나이·신기능·간기능 관련 위험 점수
abx_risk = {
    "Tazoferan(R) 4.5g": {"age": 3, "creatinine": 4, "hepatic": 2},
    "cefaZOLin 1g": {"age": 2, "creatinine": 3, "hepatic": 1},
    "Azithromycin 250mg": {"age": 1, "creatinine": 0, "hepatic": 3},
    "cefTRIAXone sod 2g": {"age": 2, "creatinine": 1, "hepatic": 3},
    "cefePIMe 1g": {"age": 3, "creatinine": 5, "hepatic": 1},
    "Amoxclan duo(R) 437.5mg/62.5mg": {"age": 1, "creatinine": 3, "hepatic": 4},
    "Meropenem 500mg": {"age": 2, "creatinine": 4, "hepatic": 2}
}

# 병원 내 주요 감염균 리스트
infectious_agents = [
    "Escherichia coli",
    "Staphylococcus aureus",
    "NTM(Nontuberculous Mycobacteria)",
    "Klebsiella pneumoniae",
    "Pseudomonas aeruginosa",
    "Enterococcus faecalis",
    "Enterococcus faecium"
]

# 그람 상태 매핑
agent_to_gram = {
    "Escherichia coli": "gram_negative",
    "Staphylococcus aureus": "gram_positive",
    "NTM(Nontuberculous Mycobacteria)": "gram_positive",
    "Klebsiella pneumoniae": "gram_negative",
    "Pseudomonas aeruginosa": "gram_negative",
    "Enterococcus faecalis": "gram_positive",
    "Enterococcus faecium": "gram_positive"
}

def validate_category(value, field_name):
    try:
        val = int(value)
        if CATEGORY_MIN <= val <= CATEGORY_MAX:
            return val
        raise ValueError(f"{field_name}는 {CATEGORY_MIN}~{CATEGORY_MAX} 사이여야 합니다.")
    except ValueError:
        raise ValueError(f"{field_name}는 숫자여야 하며 {CATEGORY_MIN}~{CATEGORY_MAX} 사이여야 합니다.")

def validate_antibiotics(input_list, valid_list, field_name):
    valid_items = [item.strip() for item in input_list if item.strip() in valid_list]
    if len(valid_items) < len([item for item in input_list if item.strip()]):
        print(f"경고: 일부 {field_name}가 유효하지 않습니다. 유효한 항목만 포함됩니다.")
    return valid_items

def get_status_score(age_category, creat_category, hepatic_category):
    return (
        validate_category(age_category, "나이 카테고리"),
        validate_category(creat_category, "신장수치 카테고리"),
        validate_category(hepatic_category, "빌리루빈 카테고리")
    )

def get_patient_states(age_score, creat_score, hepatic_score):
    states = []
    if age_score >= 4:
        states.append("extreme_age")
    if creat_score >= 4:
        states.append("extreme_creatinine")
    if hepatic_score >= 4:
        states.append("extreme_hepatic")
    return states

def recommend_antibiotics(patient):
    gram = patient['gram_status']
    agent = patient['infectious_agent']
    allergy = set(patient.get('allergy', []))
    drug_inter = set(patient.get('drug_interactions', []))
    suscept = patient.get('susceptibility', {}).get(agent, {})
    log = ["━━━━━━━━━━━━━━━━━━━━━━━", f"🔍 환자 상태 요약: {agent} ({gram})"]

    # 환자 상태
    age_score, creat_score, hepatic_score = get_status_score(
        patient['age_category'], 
        patient['renal_function']['creatinine_category'], 
        patient['hepatic_function']['bilirubin_category']
    )
    patient_states = get_patient_states(age_score, creat_score, hepatic_score)
    state_str = ", ".join(patient_states) if patient_states else "없음"
    log.append(f"  - 극도의 상태: {state_str}")
    log.append(f"  - 알러지: {', '.join(allergy) if allergy else '없음'}")
    log.append(f"  - 약물 상호작용: {', '.join(drug_inter) if drug_inter else '없음'}")
    log.append(f"  - 감수성: {', '.join([f'{k}: {v}' for k, v in suscept.items()]) if suscept else '없음'}")

    # 1. Gram 상태 필터링
    candidates = [abx for abx in abx_nodes if gram in abx_to_gram[abx]]
    eliminated_gram = [f"  · {abx}: {gram}에 효과 없음" for abx in abx_nodes if abx not in candidates]
    log.append("🔹 1. Gram 상태 필터링")
    log.append(f"  - 후보: {', '.join(candidates) if candidates else '없음'}")
    log.append("  ⮩ [제외]: " + (", ".join(eliminated_gram) if eliminated_gram else "없음"))

    # 2. 독성 필터링
    filtered = []
    eliminated_state = []
    for abx in candidates:
        toxic_reasons = []
        if "extreme_age" in patient_states and abx_risk[abx]['age'] >= 4:
            toxic_reasons.append("extreme_age")
        if "extreme_creatinine" in patient_states and abx_risk[abx]['creatinine'] >= 4:
            toxic_reasons.append("extreme_creatinine")
        if "extreme_hepatic" in patient_states and abx_risk[abx]['hepatic'] >= 4:
            toxic_reasons.append("extreme_hepatic")
        if toxic_reasons:
            eliminated_state.append(f"  · {abx}: 독성 위험 ({', '.join(toxic_reasons)})")
        else:
            filtered.append(abx)
    log.append("🔹 2. 독성 필터링")
    log.append(f"  - 후보: {', '.join(filtered) if filtered else '없음'}")
    log.append("  ⮩ [제외]: " + (", ".join(eliminated_state) if eliminated_state else "없음"))

    # 3. 알러지/약물 상호작용 필터링
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
    log.append(f"  - 후보: {', '.join(filtered2) if filtered2 else '없음'}")
    log.append("  ⮩ [제외]: " + (", ".join(eliminated_allergy) if eliminated_allergy else "없음"))

    # 4. 감수성 필터링
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

    # 5. TOPSIS 순위 계산
    if final_candidates:
        log.append("🔹 5. TOPSIS 부작용 순위")
        data = np.array([
            [get_onehot(age_score * abx_risk[abx]['age']),
             get_onehot(creat_score * abx_risk[abx]['creatinine']),
             get_onehot(hepatic_score * abx_risk[abx]['hepatic'])]
            for abx in final_candidates
        ])
        ideal, anti_ideal = np.array([0, 0, 0]), np.array([MAX_ONEHOT, MAX_ONEHOT, MAX_ONEHOT])
        dist_to_ideal = np.linalg.norm(data - ideal, axis=1)
        dist_to_anti = np.linalg.norm(data - anti_ideal, axis=1)
        Ci = dist_to_anti / (dist_to_ideal + dist_to_anti + 1e-9)
        sorted_idx = np.argsort(-Ci)
        topsis_result = [f"{final_candidates[i]} (Ci={Ci[i]:.3f})" for i in sorted_idx]
        log.append("  - 순위:")
        log.extend([f"    · {res}" for res in topsis_result])
    else:
        log.append("⚠ 추천 항생제 없음")
        topsis_result = []

    log.append("━━━━━━━━━━━━━━━━━━━━━━━")
    return final_candidates, log

def recommend_empiric_antibiotics(patient):
    log = ["━━━━━━━━━━━━━━━━━━━━━━━", "🔍 Empiric 항생제 추천 (감염균 미확인)"]
    age_score, creat_score, hepatic_score = get_status_score(
        patient['age_category'], 
        patient['renal_function']['creatinine_category'], 
        patient['hepatic_function']['bilirubin_category']
    )
    patient_states = get_patient_states(age_score, creat_score, hepatic_score)
    allergy = set(patient.get('allergy', []))
    drug_inter = set(patient.get('drug_interactions', []))
    state_str = ", ".join(patient_states) if patient_states else "없음"
    log.append(f"  - 극도의 상태: {state_str}")
    log.append(f"  - 알러지: {', '.join(allergy) if allergy else '없음'}")
    log.append(f"  - 약물 상호작용: {', '.join(drug_inter) if drug_inter else '없음'}")

    # 1. 초기 후보 (광범위 항생제 우선 고려)
    candidates = [abx for abx in abx_nodes if len(abx_to_gram[abx]) == 2]  # Broad-spectrum 우선
    candidates += [abx for abx in abx_nodes if abx not in candidates]  # 나머지 추가
    log.append("🔹 1. 초기 후보")
    log.append(f"  - 후보: {', '.join(candidates)}")

    # 2. 독성 필터링
    filtered = []
    eliminated_state = []
    for abx in candidates:
        toxic_reasons = []
        if "extreme_age" in patient_states and abx_risk[abx]['age'] >= 4:
            toxic_reasons.append("extreme_age")
        if "extreme_creatinine" in patient_states and abx_risk[abx]['creatinine'] >= 4:
            toxic_reasons.append("extreme_creatinine")
        if "extreme_hepatic" in patient_states and abx_risk[abx]['hepatic'] >= 4:
            toxic_reasons.append("extreme_hepatic")
        if toxic_reasons:
            eliminated_state.append(f"  · {abx}: 독성 위험 ({', '.join(toxic_reasons)})")
        else:
            filtered.append(abx)
    log.append("🔹 2. 독성 필터링")
    log.append(f"  - 후보: {', '.join(filtered) if filtered else '없음'}")
    log.append("  ⮩ [제외]: " + (", ".join(eliminated_state) if eliminated_state else "없음"))

    # 3. 알러지/약물 상호작용 필터링
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

    # 4. TOPSIS 순위 계산
    if final_candidates:
        log.append("🔹 4. TOPSIS 부작용 순위")
        data = np.array([
            [get_onehot(age_score * abx_risk[abx]['age']),
             get_onehot(creat_score * abx_risk[abx]['creatinine']),
             get_onehot(hepatic_score * abx_risk[abx]['hepatic'])]
            for abx in final_candidates
        ])
        ideal, anti_ideal = np.array([0, 0, 0]), np.array([MAX_ONEHOT, MAX_ONEHOT, MAX_ONEHOT])
        dist_to_ideal = np.linalg.norm(data - ideal, axis=1)
        dist_to_anti = np.linalg.norm(data - anti_ideal, axis=1)
        Ci = dist_to_anti / (dist_to_ideal + dist_to_anti + 1e-9)
        sorted_idx = np.argsort(-Ci)
        topsis_result = [f"{final_candidates[i]} (Ci={Ci[i]:.3f})" for i in sorted_idx]
        log.append("  - 순위:")
        log.extend([f"    · {res}" for res in topsis_result])
    else:
        log.append("⚠ 추천 항생제 없음")
        topsis_result = []

    log.append("━━━━━━━━━━━━━━━━━━━━━━━")
    return final_candidates, log

def get_user_input():
    print("=== 환자 데이터 입력 ===")
    patient_id = input("환자 ID를 입력하세요: ").strip()

    # 감염균 확인 여부
    pathogen_known = input("감염균을 아나요? (y/n): ").strip().lower() == 'y'

    # 공통 입력
    print("나이 카테고리 (1~5): 1: 10세 이하, 2: 11~30세, 3: 31~50세, 4: 51~60세, 5: 61세 이상")
    age_category = validate_category(input("선택: "), "나이 카테고리")

    print("\n신장수치(크레아티닌) 카테고리 (1~5): 1: 0.9 이하, 2: 0.9~1.2, 3: 1.2~1.8, 4: 1.8~2.3, 5: 2.3 초과")
    creatinine_category = validate_category(input("선택: "), "신장수치 카테고리")

    print("\n빌리루빈 카테고리 (1~5): 1: 0.9 이하, 2: 0.9~1.2, 3: 1.2~1.8, 4: 1.8~2.3, 5: 2.3 초과")
    bilirubin_category = validate_category(input("선택: "), "빌리루빈 카테고리")

    print("\n알러지 항생제 (쉼표로 구분, 없으면 엔터, 유효 항생제: " + ", ".join(abx_nodes) + "):")
    allergy_input = input().strip()
    allergy = validate_antibiotics([a.strip() for a in allergy_input.split(",") if a.strip()], abx_nodes, "알러지")

    print("\n약물 상호작용 항생제 (쉼표로 구분, 없으면 엔터, 유효 항생제: " + ", ".join(abx_nodes) + "):")
    drug_inter_input = input().strip()
    drug_inter = validate_antibiotics([d.strip() for d in drug_inter_input.split(",") if d.strip()], abx_nodes, "약물 상호작용")

    patient = {
        "patient_id": patient_id,
        "age_category": age_category,
        "renal_function": {"creatinine_category": creatinine_category},
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
        susceptibility = {}
        print("\n항생제별 감수성(S, I, R, 기본값 S):")
        for abx in abx_nodes:
            sir = input(f"{abx}: ").strip().upper()
            susceptibility[abx] = "S" if sir == "" else sir
            if sir not in ["", "S", "I", "R"]:
                print(f"경고: {abx}의 감수성 {sir}은 유효하지 않음. 기본값 S로 설정.")
        patient["susceptibility"] = {infectious_agent: susceptibility}

    return patient

def main():
    print("=== 항생제 추천 시스템 ===")
    while True:
        patient = get_user_input()
        print("\n=== 환자 정보 요약 ===")
        print(f"ID: {patient['patient_id']}")
        print(f"나이 카테고리: {patient['age_category']}")
        print(f"신장수치 카테고리: {patient['renal_function']['creatinine_category']}")
        print(f"빌리루빈 카테고리: {patient['hepatic_function']['bilirubin_category']}")
        print(f"알러지: {', '.join(patient['allergy']) if patient['allergy'] else '없음'}")
        print(f"약물 상호작용: {', '.join(patient['drug_interactions']) if patient['drug_interactions'] else '없음'}")
        if patient['pathogen_known']:
            print(f"감염균: {patient['infectious_agent']} ({patient['gram_status']})")
            print("감수성:")
            for abx, sir in patient['susceptibility'][patient['infectious_agent']].items():
                print(f"  {abx}: {sir}")
        else:
            print("감염균: 미확인")

        print("\n=== 추천 결과 ===")
        if patient['pathogen_known']:
            result, log = recommend_antibiotics(patient)
        else:
            result, log = recommend_empir
