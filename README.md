#항생제추천시스템

# streamlit_app.py
항생제 리스트 불러오고 (주로 쓰인 항생제 CDW 에서 가져옴)
항생제별 gram +/- 로 effective 를 설정 (바꿀 필요 잇음)
A/B 별 risk score 를 임의로 결정
Patient 별 risk score 를 range 에 따라 결정
A/B + Patient 별 risk score 를 combination 하여 결정
항생제 추천 

# Patient Example.json
예시 가짜 환자 json 형식 데이터
비정형 데이터 기반:
    "susceptibility": {
      "Escherichia coli": {
        "Tazoferan(R) 4.5g": "S",
        "cefaZOLin 1g": "R",
        "Azithromycin 250mg": "R",
        "cefTRIAXone sod 2g": "S",
        "cefePIMe 1g": "S",
        "Amoxclan duo(R) 437.5mg/62.5mg": "I",
        "Meropenem 500mg": "S"
- 환자의 항생제 감수성 정보 추출
