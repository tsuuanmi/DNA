//! Test support: the output contract of `dna-eval run`, `dna.eval_run/v1`
//! (`docs/reference/schemas/eval-run-v1.schema.json`, PROP-0011), read so the
//! tests check what the harness writes against it. A table's comment lines,
//! header, cell count, and enumerated cells are checked here; the cell patterns
//! and the whole directory are checked by `validate_result_schemas.py`.

use serde_json::{Map, Value};

/// The schema, as committed.
const SCHEMA: &str = include_str!("../../../docs/reference/schemas/eval-run-v1.schema.json");

/// A contract check's outcome.
pub(crate) type Checked = Result<(), Box<dyn std::error::Error>>;

fn schema() -> Result<Value, Box<dyn std::error::Error>> {
    Ok(serde_json::from_str(SCHEMA)?)
}

/// `node`, or the definition its local `$ref` names.
fn resolve<'a>(schema: &'a Value, node: &'a Value) -> Result<&'a Value, String> {
    match node.get("$ref").and_then(Value::as_str) {
        Some(reference) => {
            let name = reference
                .strip_prefix("#/$defs/")
                .ok_or_else(|| format!("unsupported reference {reference}"))?;
            schema["$defs"]
                .get(name)
                .ok_or_else(|| format!("no definition {name}"))
        }
        None => Ok(node),
    }
}

/// The definition of the file `name`.
fn file<'a>(schema: &'a Value, name: &str) -> Result<&'a Value, String> {
    let node = schema["properties"]
        .get(name)
        .ok_or_else(|| format!("the contract names no file {name}"))?;
    resolve(schema, node)
}

/// Checks the table the harness writes as `name`: its comment lines, its
/// header, and the cells of every row.
pub(crate) fn check_table(name: &str, text: &str) -> Checked {
    let schema = schema()?;
    let table = file(&schema, name)?;
    let properties = &table["properties"];
    let mut lines = text
        .strip_suffix('\n')
        .ok_or_else(|| format!("{name} does not end in a newline"))?
        .split('\n')
        .peekable();
    let mut comments = 0;
    while lines.next_if(|line| line.starts_with("# ")).is_some() {
        comments += 1;
    }
    let comment_schema = resolve(&schema, &properties["comments"])?;
    let bounds = |key: &str| comment_schema.get(key).and_then(Value::as_u64);
    if bounds("minItems").is_some_and(|minimum| comments < minimum)
        || bounds("maxItems").is_some_and(|maximum| comments > maximum)
    {
        return Err(format!("{name} has {comments} comment lines").into());
    }
    let header = lines
        .next()
        .ok_or_else(|| format!("{name} has no header"))?
        .split('\t')
        .collect::<Vec<_>>();
    let columns = &properties["columns"];
    let fixed = match columns.get("const") {
        Some(Value::Array(names)) => names.iter().collect::<Vec<_>>(),
        _ => columns["prefixItems"]
            .as_array()
            .map(|items| items.iter().map(|item| &item["const"]).collect())
            .unwrap_or_default(),
    };
    let fixed = fixed
        .iter()
        .map(|name| name.as_str().unwrap_or_default())
        .collect::<Vec<_>>();
    let open = columns.get("items").is_some();
    if header.get(..fixed.len()) != Some(&fixed[..]) || (!open && header.len() != fixed.len()) {
        return Err(format!("{name} header {header:?} is not the contract's {fixed:?}").into());
    }
    let row = &properties["rows"]["items"];
    for (number, line) in lines.enumerate() {
        let cells = line.split('\t').collect::<Vec<_>>();
        if cells.len() != header.len() {
            return Err(format!("{name} row {} has {} cells", number + 1, cells.len()).into());
        }
        for (column, cell) in header.iter().zip(&cells) {
            let cell_schema = row["properties"]
                .get(*column)
                .unwrap_or(&row["additionalProperties"]);
            check_value(
                &schema,
                cell_schema,
                &Value::from(*cell),
                &format!("{name}:{column}"),
            )?;
        }
    }
    Ok(())
}

/// Checks files the harness writes, by name: every table, and that the
/// contract names each file.
pub(crate) fn check_files<'a>(files: impl IntoIterator<Item = (&'a str, &'a str)>) -> Checked {
    let schema = schema()?;
    for (name, text) in files {
        file(&schema, name)?;
        if std::path::Path::new(name)
            .extension()
            .is_some_and(|extension| extension == "tsv")
        {
            check_table(name, text)?;
        }
    }
    Ok(())
}

/// Checks the JSON file the harness writes as `name`.
pub(crate) fn check_document(name: &str, value: &Value) -> Checked {
    let schema = schema()?;
    let node = file(&schema, name)?;
    check_value(&schema, node, value, name)
}

/// Checks `value` against `node`: its type, `const`, `enum`, bounds, and an
/// object's fields, recursively; patterns are not checked.
fn check_value(schema: &Value, node: &Value, value: &Value, at: &str) -> Checked {
    let node = resolve(schema, node)?;
    let fail = |what: &str| -> Checked { Err(format!("{at}: {what}: {value}").into()) };
    if node.get("const").is_some_and(|expected| expected != value) {
        return fail("not the constant");
    }
    if let Some(Value::Array(allowed)) = node.get("enum")
        && !allowed.contains(value)
    {
        return fail("not an allowed value");
    }
    if let Some(types) = node.get("type") {
        let types = match types {
            Value::Array(types) => types.iter().filter_map(Value::as_str).collect::<Vec<_>>(),
            other => other.as_str().into_iter().collect(),
        };
        let found = match value {
            Value::Null => "null",
            Value::Bool(_) => "boolean",
            Value::Number(number) if number.is_i64() || number.is_u64() => "integer",
            Value::Number(_) => "number",
            Value::String(_) => "string",
            Value::Array(_) => "array",
            Value::Object(_) => "object",
        };
        // An integer is also a number.
        let number = found == "integer" && types.contains(&"number");
        if !number && !types.contains(&found) {
            return fail("of another type");
        }
    }
    let bound = |key: &str| node.get(key).and_then(Value::as_f64);
    if let Some(number) = value.as_f64()
        && (bound("minimum").is_some_and(|minimum| number < minimum)
            || bound("maximum").is_some_and(|maximum| number > maximum))
    {
        return fail("out of bounds");
    }
    match value {
        Value::Object(fields) => check_fields(schema, node, fields, at),
        Value::Array(items) => match node.get("items") {
            Some(item) if item.is_object() => {
                items.iter().enumerate().try_for_each(|(index, value)| {
                    check_value(schema, item, value, &format!("{at}[{index}]"))
                })
            }
            _ => Ok(()),
        },
        _ => Ok(()),
    }
}

/// Checks an object's required and allowed fields, and each field's value.
fn check_fields(schema: &Value, node: &Value, fields: &Map<String, Value>, at: &str) -> Checked {
    if let Some(Value::Array(required)) = node.get("required")
        && let Some(missing) = required
            .iter()
            .filter_map(Value::as_str)
            .find(|name| !fields.contains_key(*name))
    {
        return Err(format!("{at}: lacks {missing}").into());
    }
    for (name, value) in fields {
        let at = format!("{at}.{name}");
        match (
            node["properties"].get(name),
            node.get("additionalProperties"),
        ) {
            (Some(field), _) | (None, Some(field @ Value::Object(_))) => {
                check_value(schema, field, value, &at)?;
            }
            (None, Some(Value::Bool(false))) => {
                return Err(format!("{at}: not in the contract").into());
            }
            (None, _) => {}
        }
    }
    Ok(())
}

/// The checks reject what the contract does not describe.
#[test]
fn rejects_what_the_contract_does_not_describe() -> Checked {
    let header = "sample_id\tfirst\tlast\tlabels\n";
    check_table("checks.tsv", header)?;
    check_table("checks.tsv", &format!("{header}S1\t1\t2\tone_read\n"))?;
    let windows = "region\tstart\tend\tsamples\tchecked\tunreviewed\tunreviewed_share\t\
                   unreviewed_upper95\n";
    check_table(
        "windows.tsv",
        &format!("{windows}HV1\t1\t10\t0\t0\t0\t\t1.0000\n"),
    )?;
    for (name, text) in [
        ("notes.tsv", header.to_owned()),
        ("checks.tsv", "sample_id\tfirst\tlast\n".to_owned()),
        ("checks.tsv", format!("{header}S1\t1\t2\n")),
        ("checks.tsv", format!("# note\n{header}")),
        ("checks.tsv", header.trim_end().to_owned()),
        (
            "windows.tsv",
            format!("{windows}HV4\t1\t10\t0\t0\t0\t\t1.0000\n"),
        ),
    ] {
        assert!(check_table(name, &text).is_err(), "{name}: {text:?}");
    }
    let mut example: Value = serde_json::from_str(include_str!(
        "../../../docs/reference/examples/eval-run-v1.example.json"
    ))?;
    let run = &mut example["run.json"];
    check_document("run.json", run)?;
    for (field, value) in [("threads", Value::from(0)), ("extra", Value::from(1))] {
        let mut changed = run.clone();
        changed[field] = value;
        assert!(check_document("run.json", &changed).is_err(), "{field}");
    }
    run["truth"]["state"] = Value::from("stale");
    assert!(check_document("run.json", run).is_err());
    Ok(())
}
