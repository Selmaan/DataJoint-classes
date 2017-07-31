function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'v1inf', 'v1inf');
end
obj = schemaObject;
end
