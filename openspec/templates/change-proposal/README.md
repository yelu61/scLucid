# OpenSpec Change Proposal Template (scLucid)

Use this template set to scaffold a new change proposal quickly.

## 1. Choose a change id
Use a unique, verb-led kebab-case id, for example:
- `add-tumor-clone-quality-metrics`
- `update-qc-adaptive-thresholding`
- `refactor-analysis-import-stability`

## 2. Create the change directory
```bash
CHANGE_ID="add-your-change-name"
mkdir -p "openspec/changes/${CHANGE_ID}/specs/<capability>"
```

## 3. Copy templates
```bash
cp openspec/templates/change-proposal/proposal.md "openspec/changes/${CHANGE_ID}/proposal.md"
cp openspec/templates/change-proposal/tasks.md "openspec/changes/${CHANGE_ID}/tasks.md"
cp openspec/templates/change-proposal/spec-delta.md "openspec/changes/${CHANGE_ID}/specs/<capability>/spec.md"
```

Optional design doc:
```bash
cp openspec/templates/change-proposal/design.md "openspec/changes/${CHANGE_ID}/design.md"
```

## 4. Fill in all placeholders
Replace `<...>` tokens with concrete scLucid details.

## 5. Validate
```bash
openspec validate "${CHANGE_ID}" --strict
```
