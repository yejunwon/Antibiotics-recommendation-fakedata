import streamlit as st
import json
import networkx as nx
import numpy as np

# One-hot encoding 매핑 (0~25 -> 0~14)
score2onehot = {
    0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 8: 7, 9: 8,
    10: 9, 12: 10, 15: 11, 16: 12, 20: 13, 25: 14
}

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

def get_status_score(age_category, creat_category, hepatic_category):
    if not (1 <= age_category <= 5):
        raise ValueError("나이 카테고리는 1~5 사이여야 합니다.")
    if not (1 <= creat_category <= 5):
        raise ValueError("신장수치 카테고리는 1~5 사이여야 합니다.")
    if not (1 <= hepatic_category <= 5):
        raise ValueError("빌리루빈 카테고리는 1~5 사이여야 합니다.")
    return age_category, creat_category, hepatic_category

# 지식 그래프 생성
gram_nodes = ['gram_positive', 'gram_negative']
state_nodes = ['extreme_age', 'extreme_creatinine', 'extreme_hepatic']
KG = nx.DiGraph()
KG.add_nodes_from(gram_nodes, bipartite=0)
KG.add_nodes_from(state_nodes, bipartite=2)
KG.add_nodes_from(abx_nodes, bipartite=1)
for abx, grams in abx_to_gram.items():
    for gram in grams:
        KG.add_edge(abx, gram, relation='effective_against')
for abx, risk in abx_risk.items():
    if risk['age'] >= 4:
        KG.add_edge(abx, 'extreme_age', relation='is_toxic_to')
    if risk['creatinine'] >= 4:
        KG.add_edge(abx, 'extreme_creatinine', relation='is_toxic_to')
    if risk['hepatic'] >= 4:
        KG.add_edge(abx, 'extreme_hepatic', relation='is_toxic_to')

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
    log = []

    # 0. 환자 상태 요약
    patient_states = get_patient_states(*get_status_score(
        patient['age_category'], 
        patient['renal_function']['creatinine_category'], 
        patient['hepatic_function']['bilirubin_category']
    ))
    state_str = ", ".join(patient_states) if patient_states else "없음"
    log.append("━━━━━━━━━━━━━━━━━━━━━━━")
    log.append(f"🔍 환자 상태 요약")
    log.append(f"  - 감염균: {agent} ({gram})")
    log.append(f"  - 극도의 상태: {state_str}")
    log.append(f"  - 알러지: {', '.join(allergy) if allergy else '없음'}")
    log.append(f"  - 약물 상호작용: {', '.join(drug_inter) if drug_inter else '없음'}")
    log.append(f"  - 감수성 데이터: {', '.join([f'{k}: {v}' for k, v in suscept.items()]) if suscept else '없음'}")
    log.append("━━━━━━━━━━━━━━━━━━━━━━━\n")

    # 1단계: Gram 상태 기반 후보
    candidates = []
    eliminated_gram = []
    for abx in abx_nodes:
        if not KG.has_edge(abx, gram):
            eliminated_gram.append(f"  · {abx}: {gram}에 효과 없음")
            continue
        candidates.append(abx)
    log.append("🔹 1단계: Gram 상태 기반 필터링")
    log.append(f"  - 고려 대상: {', '.join(candidates) if candidates else '없음'}")
    if eliminated_gram:
        log.append("  ⮩ [제외 항목]")
        log.extend(eliminated_gram)
    else:
        log.append("  ⮩ [제외 항목]: 없음")
    log.append("")

    # 2단계: 환자 상태(나이, 신장, 간기능) 기반 독성 필터링
    filtered = []
    eliminated_state = []
    for abx in candidates:
        toxic = False
        toxic_reasons = []
        for state in patient_states:
            if KG.has_edge(abx, state) and KG[abx][state].get('relation') == 'is_toxic_to':
                toxic = True
                toxic_reasons.append(state)
        if toxic:
            eliminated_state.append(f"  · {abx}: 독성 위험 ({', '.join(toxic_reasons)})")
        else:
            filtered.append(abx)
    log.append("🔹 2단계: 환자 상태 기반 독성 필터링")
    log.append(f"  - 고려 대상: {', '.join(filtered) if filtered else '없음'}")
    if eliminated_state:
        log.append("  ⮩ [제외 항목]")
        log.extend(eliminated_state)
    else:
        log.append("  ⮩ [제외 항목]: 없음")
    log.append("")

    # 3단계: 알러지/약물 상호작용 필터링
    filtered2 = []
    eliminated_allergy_inter = []
    reason_dict = {}
    for abx in filtered:
        reasons = []
        if abx in allergy:
            reasons.append("알러지")
        if abx in drug_inter:
            reasons.append("약물 상호작용")
        if reasons:
            eliminated_allergy_inter.append(f"  · {abx}: {', '.join(reasons)}")
            reason_dict[abx] = ', '.join(reasons)
        else:
            filtered2.append(abx)
    log.append("🔹 3단계: 알러지 및 약물 상호작용 필터링")
    log.append(f"  - 고려 대상: {', '.join(filtered2) if filtered2 else '없음'}")
    if eliminated_allergy_inter:
        log.append("  ⮩ [제외 항목]")
        log.extend(eliminated_allergy_inter)
    else:
        log.append("  ⮩ [제외 항목]: 없음")
    log.append("")

    # 4단계: 감수성(S) 필터링, 없으면 I 고려
    final_candidates = []
    eliminated_suscept = []
    for abx in filtered2:
        if suscept.get(abx, None) == "S":
            final_candidates.append(abx)
        else:
            eliminated_suscept.append(f"  · {abx}: 감수성 미달 (S 아님, 현재: {suscept.get(abx, '없음')})")
    if not final_candidates:
        for abx in filtered2:
            if suscept.get(abx, None) == "I":
                final_candidates.append(abx)
                log.append(f"  · {abx}: 감수성 I (S 없으므로 포함)")
            else:
                eliminated_suscept.append(f"  · {abx}: 감수성 미달 (S/I 아님, 현재: {suscept.get(abx, '없음')})")
    log.append("🔹 4단계: 감수성(S 또는 I) 필터링")
    log.append(f"  - 최종 후보: {', '.join(final_candidates) if final_candidates else '없음'}")
    if eliminated_suscept:
        log.append("  ⮩ [제외 항목]")
        log.extend(eliminated_suscept)
    else:
        log.append("  ⮩ [제외 항목]: 없음")
    log.append("")

    # 5단계: 최종 추천
    log.append("🔹 5단계: 최종 추천 항생제")
    if final_candidates:
        for abx in final_candidates:
            log.append(f"  · {abx}")
    else:
        log.append("  ⚠ 추천 항생제 없음")
    log.append("")

    # Toxicity Score 계산 및 로그
    age_score, creat_score, hepatic_score = get_status_score(
        patient['age_category'], 
        patient['renal_function']['creatinine_category'], 
        patient['hepatic_function']['bilirubin_category']
    )
    tox_info = ["━━━━━━━━━━━━━━━━━━━━━━━", "💊 항생제 Toxicity Score"]
    tox_info.append(f"  - 환자 점수: 나이={age_score}, 신장={creat_score}, 간기능={hepatic_score}")
    for abx in final_candidates:
        a_risk = abx_risk[abx]['age']
        c_risk = abx_risk[abx]['creatinine']
        h_risk = abx_risk[abx]['hepatic']
        age_raw = age_score * a_risk
        creat_raw = creat_score * c_risk
        hepatic_raw = hepatic_score * h_risk
        age_tox = get_onehot(age_raw)
        creat_tox = get_onehot(creat_raw)
        hepatic_tox = get_onehot(hepatic_raw)
        tox_info.append(f"  · {abx:12}:")
        tox_info.append(f"     - 나이: {age_score} × {a_risk} = {age_raw} → one-hot: {age_tox}")
        tox_info.append(f"     - 신장: {creat_score} × {c_risk} = {creat_raw} → one-hot: {creat_tox}")
        tox_info.append(f"     - 간기능: {hepatic_score} × {h_risk} = {hepatic_raw} → one-hot: {hepatic_tox}")
    log += tox_info + [""]

    # TOPSIS 순위 계산 및 로그
    if final_candidates:
        log.append("🔹 TOPSIS 부작용 순위 계산")
        A_list, C_list, H_list = [], [], []
        for abx in final_candidates:
            a_risk = abx_risk[abx]['age']
            c_risk = abx_risk[abx]['creatinine']
            h_risk = abx_risk[abx]['hepatic']
            A_list.append(get_onehot(age_score * a_risk))
            C_list.append(get_onehot(creat_score * c_risk))
            H_list.append(get_onehot(hepatic_score * h_risk))
        data = np.array(list(zip(A_list, C_list, H_list)))
        log.append("  - 점수 배열:")
        for i, abx in enumerate(final_candidates):
            log.append(f"    · {abx}: [나이={data[i,0]}, 신장={data[i,1]}, 간기능={data[i,2]}]")
        ideal = np.array([0, 0, 0])
        anti_ideal = np.array([14, 14, 14])
        log.append(f"  - 이상점: [0, 0, 0]")
        log.append(f"  - 반이상점: [14, 14, 14]")
        dist_to_ideal = np.linalg.norm(data - ideal, axis=1)
        dist_to_anti = np.linalg.norm(data - anti_ideal, axis=1)
        log.append("  - 이상점까지 거리:")
        for i, abx in enumerate(final_candidates):
            log.append(f"    · {abx}: {dist_to_ideal[i]:.3f}")
        log.append("  - 반이상점까지 거리:")
        for i, abx in enumerate(final_candidates):
            log.append(f"    · {abx}: {dist_to_anti[i]:.3f}")
        Ci = dist_to_anti / (dist_to_ideal + dist_to_anti + 1e-9)
        if np.any(np.isnan(Ci)):
            log.append("  ⚠ TOPSIS 계산 중 오류: 유효하지 않은 Ci 값")
            topsis_result = final_candidates
        else:
            sorted_idx = np.argsort(-Ci)
            topsis_result = [f"{final_candidates[i]} (Ci={Ci[i]:.3f})" for i in sorted_idx]
            log.append("  - Ci 값 (높을수록 부작용 적음):")
            for i, abx in enumerate(final_candidates):
                log.append(f"    · {abx}: {Ci[i]:.3f}")
            log.append("  - 최종 순위:")
            for rec in topsis_result:
                log.append(f"    · {rec}")
    else:
        topsis_result = []
        log.append("🔹 TOPSIS 부작용 순위 계산")
        log.append("  ⚠ 추천 항생제 없음, TOPSIS 계산 생략")
    log.append("━━━━━━━━━━━━━━━━━━━━━━━")
    return final_candidates, log

def main():
    st.title("항생제 추천 시스템")

    # 환자 데이터 입력
    st.header("환자 데이터 입력")
    patient_id = st.text_input("환자 ID", value="환자1")

    # 나이 카테고리
    age_options = {
        "10세 이하": 1,
        "11~30세": 2,
        "31~50세": 3,
        "51~60세": 4,
        "61세 이상": 5
    }
    age_category_label = st.selectbox("나이 카테고리", list(age_options.keys()))
    age_category = age_options[age_category_label]

    # 신장수치 카테고리
    creatinine_options = {
        "0.9 이하": 1,
        "0.9~1.2": 2,
        "1.2~1.8": 3,
        "1.8~2.3": 4,
        "2.3 초과": 5
    }
    creatinine_category_label = st.selectbox("신장수치(크레아티닌) 카테고리", list(creatinine_options.keys()))
    creatinine_category = creatinine_options[creatinine_category_label]

    # 빌리루빈 카테고리
    bilirubin_options = {
        "0.9 이하": 1,
        "0.9~1.2": 2,
        "1.2~1.8": 3,
        "1.8~2.3": 4,
        "2.3 초과": 5
    }
    bilirubin_category_label = st.selectbox("빌리루빈 카테고리", list(bilirubin_options.keys()))
    bilirubin_category = bilirubin_options[bilirubin_category_label]

    # 감염균 선택
    infectious_agent = st.selectbox("감염균", infectious_agents)
    gram_status = agent_to_gram[infectious_agent]

    # 알러지 선택
    allergy = st.multiselect("알러지 항생제 (복수 선택 가능)", abx_nodes, default=[])

    # 약물 상호작용 선택
    drug_inter = st.multiselect("약물 상호작용 항생제 (복수 선택 가능)", abx_nodes, default=[])

    # 감수성 입력
    st.subheader("항생제별 감수성 (S, I, R)")
    susceptibility = {}
    for abx in abx_nodes:
        sir = st.selectbox(f"{abx} 감수성", ["S", "I", "R"], index=0, key=f"suscept_{abx}")
        susceptibility[abx] = sir

    # 환자 데이터 구성
    patient = {
        "patient_id": patient_id,
        "age_category": age_category,
        "renal_function": {"creatinine_category": creatinine_category},
        "hepatic_function": {"bilirubin_category": bilirubin_category},
        "infectious_agent": infectious_agent,
        "gram_status": gram_status,
        "allergy": allergy,
        "drug_interactions": drug_inter,
        "susceptibility": {infectious_agent: susceptibility}
    }

    # 환자 정보 요약
    st.header("환자 정보 요약")
    st.write(f"**ID**: {patient['patient_id']}")
    st.write(f"**나이 카테고리**: {age_category_label} ({age_category})")
    st.write(f"**신장수치 카테고리**: {creatinine_category_label} ({creatinine_category})")
    st.write(f"**빌리루빈 카테고리**: {bilirubin_category_label} ({bilirubin_category})")
    st.write(f"**감염균**: {patient['infectious_agent']} ({patient['gram_status']})")
    st.write(f"**알러지**: {', '.join(patient['allergy']) if patient['allergy'] else '없음'}")
    st.write(f"**약물 상호작용**: {', '.join(patient['drug_interactions']) if patient['drug_interactions'] else '없음'}")
    st.write("**감수성**:")
    for abx, sir in patient['susceptibility'][patient['infectious_agent']].items():
        st.write(f"  - {abx}: {sir}")

    # 추천 버튼
    if st.button("항생제 추천"):
        st.header("추천 결과")
        result, log = recommend_antibiotics(patient)
        if result:
            st.success("추천 항생제:")
            for abx in result:
                st.write(f"- {abx}")
        else:
            st.warning("⚠ 추천 항생제가 없습니다.")

        st.subheader("추천 Reasoning Log")
        st.text("\n".join(log))

if __name__ == "__main__":
    main()