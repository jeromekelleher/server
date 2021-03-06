# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: proto/ga4gh/references.proto

import sys
_b=sys.version_info[0]<3 and (lambda x:x) or (lambda x:x.encode('latin1'))
from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from google.protobuf import reflection as _reflection
from google.protobuf import symbol_database as _symbol_database
from google.protobuf import descriptor_pb2
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()




DESCRIPTOR = _descriptor.FileDescriptor(
  name='proto/ga4gh/references.proto',
  package='ga4gh',
  syntax='proto3',
  serialized_pb=_b('\n\x1cproto/ga4gh/references.proto\x12\x05ga4gh\"\x90\x01\n\tReference\x12\n\n\x02id\x18\x01 \x01(\t\x12\x0e\n\x06length\x18\x02 \x01(\x03\x12\x13\n\x0bmd5checksum\x18\x03 \x01(\t\x12\x0c\n\x04name\x18\x04 \x01(\t\x12\x12\n\nsource_uri\x18\x05 \x01(\t\x12\x19\n\x11source_accessions\x18\x06 \x03(\t\x12\x15\n\rncbi_taxon_id\x18\x07 \x01(\x05\"\xb6\x01\n\x0cReferenceSet\x12\n\n\x02id\x18\x01 \x01(\t\x12\x15\n\rreference_ids\x18\x02 \x03(\t\x12\x13\n\x0bmd5checksum\x18\x03 \x01(\t\x12\x15\n\rncbi_taxon_id\x18\x04 \x01(\x05\x12\x13\n\x0b\x64\x65scription\x18\x05 \x01(\t\x12\x13\n\x0b\x61ssembly_id\x18\x06 \x01(\t\x12\x12\n\nsource_uri\x18\x07 \x01(\t\x12\x19\n\x11source_accessions\x18\x08 \x03(\t\"\x82\x01\n\x1aSearchReferenceSetsRequest\x12\x14\n\x0cmd5checksums\x18\x01 \x03(\t\x12\x12\n\naccessions\x18\x02 \x03(\t\x12\x13\n\x0b\x61ssembly_id\x18\x03 \x01(\t\x12\x12\n\npage_token\x18\x04 \x01(\t\x12\x11\n\tpage_size\x18\x05 \x01(\x05\"c\n\x1bSearchReferenceSetsResponse\x12+\n\x0ereference_sets\x18\x01 \x03(\x0b\x32\x13.ga4gh.ReferenceSet\x12\x17\n\x0fnext_page_token\x18\x02 \x01(\t\"2\n\x16GetReferenceSetRequest\x12\x18\n\x10reference_set_id\x18\x01 \x01(\t\"\x84\x01\n\x17SearchReferencesRequest\x12\x14\n\x0cmd5checksums\x18\x01 \x03(\t\x12\x12\n\naccessions\x18\x02 \x03(\t\x12\x18\n\x10reference_set_id\x18\x03 \x01(\t\x12\x12\n\npage_token\x18\x04 \x01(\t\x12\x11\n\tpage_size\x18\x05 \x01(\x05\"Y\n\x18SearchReferencesResponse\x12$\n\nreferences\x18\x01 \x03(\x0b\x32\x10.ga4gh.Reference\x12\x17\n\x0fnext_page_token\x18\x02 \x01(\t\"+\n\x13GetReferenceRequest\x12\x14\n\x0creference_id\x18\x01 \x01(\t\"k\n\x10ListBasesRequest\x12\x14\n\x0creference_id\x18\x01 \x01(\t\x12\r\n\x05start\x18\x02 \x01(\x03\x12\x0b\n\x03\x65nd\x18\x03 \x01(\x03\x12\x12\n\npage_token\x18\x04 \x01(\t\x12\x11\n\tpage_size\x18\x05 \x01(\x05\"N\n\x11ListBasesResponse\x12\x0e\n\x06offset\x18\x01 \x01(\x03\x12\x10\n\x08sequence\x18\x02 \x01(\t\x12\x17\n\x0fnext_page_token\x18\x03 \x01(\t2\x8a\x03\n\x10ReferenceService\x12\\\n\x13SearchReferenceSets\x12!.ga4gh.SearchReferenceSetsRequest\x1a\".ga4gh.SearchReferenceSetsResponse\x12\x45\n\x0fGetReferenceSet\x12\x1d.ga4gh.GetReferenceSetRequest\x1a\x13.ga4gh.ReferenceSet\x12S\n\x10SearchReferences\x12\x1e.ga4gh.SearchReferencesRequest\x1a\x1f.ga4gh.SearchReferencesResponse\x12<\n\x0cGetReference\x12\x1a.ga4gh.GetReferenceRequest\x1a\x10.ga4gh.Reference\x12>\n\tListBases\x12\x17.ga4gh.ListBasesRequest\x1a\x18.ga4gh.ListBasesResponseb\x06proto3')
)
_sym_db.RegisterFileDescriptor(DESCRIPTOR)




_REFERENCE = _descriptor.Descriptor(
  name='Reference',
  full_name='ga4gh.Reference',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='id', full_name='ga4gh.Reference.id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='length', full_name='ga4gh.Reference.length', index=1,
      number=2, type=3, cpp_type=2, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='md5checksum', full_name='ga4gh.Reference.md5checksum', index=2,
      number=3, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='name', full_name='ga4gh.Reference.name', index=3,
      number=4, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='source_uri', full_name='ga4gh.Reference.source_uri', index=4,
      number=5, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='source_accessions', full_name='ga4gh.Reference.source_accessions', index=5,
      number=6, type=9, cpp_type=9, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='ncbi_taxon_id', full_name='ga4gh.Reference.ncbi_taxon_id', index=6,
      number=7, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=40,
  serialized_end=184,
)


_REFERENCESET = _descriptor.Descriptor(
  name='ReferenceSet',
  full_name='ga4gh.ReferenceSet',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='id', full_name='ga4gh.ReferenceSet.id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='reference_ids', full_name='ga4gh.ReferenceSet.reference_ids', index=1,
      number=2, type=9, cpp_type=9, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='md5checksum', full_name='ga4gh.ReferenceSet.md5checksum', index=2,
      number=3, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='ncbi_taxon_id', full_name='ga4gh.ReferenceSet.ncbi_taxon_id', index=3,
      number=4, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='description', full_name='ga4gh.ReferenceSet.description', index=4,
      number=5, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='assembly_id', full_name='ga4gh.ReferenceSet.assembly_id', index=5,
      number=6, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='source_uri', full_name='ga4gh.ReferenceSet.source_uri', index=6,
      number=7, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='source_accessions', full_name='ga4gh.ReferenceSet.source_accessions', index=7,
      number=8, type=9, cpp_type=9, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=187,
  serialized_end=369,
)


_SEARCHREFERENCESETSREQUEST = _descriptor.Descriptor(
  name='SearchReferenceSetsRequest',
  full_name='ga4gh.SearchReferenceSetsRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='md5checksums', full_name='ga4gh.SearchReferenceSetsRequest.md5checksums', index=0,
      number=1, type=9, cpp_type=9, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='accessions', full_name='ga4gh.SearchReferenceSetsRequest.accessions', index=1,
      number=2, type=9, cpp_type=9, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='assembly_id', full_name='ga4gh.SearchReferenceSetsRequest.assembly_id', index=2,
      number=3, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_token', full_name='ga4gh.SearchReferenceSetsRequest.page_token', index=3,
      number=4, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_size', full_name='ga4gh.SearchReferenceSetsRequest.page_size', index=4,
      number=5, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=372,
  serialized_end=502,
)


_SEARCHREFERENCESETSRESPONSE = _descriptor.Descriptor(
  name='SearchReferenceSetsResponse',
  full_name='ga4gh.SearchReferenceSetsResponse',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='reference_sets', full_name='ga4gh.SearchReferenceSetsResponse.reference_sets', index=0,
      number=1, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='next_page_token', full_name='ga4gh.SearchReferenceSetsResponse.next_page_token', index=1,
      number=2, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=504,
  serialized_end=603,
)


_GETREFERENCESETREQUEST = _descriptor.Descriptor(
  name='GetReferenceSetRequest',
  full_name='ga4gh.GetReferenceSetRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='reference_set_id', full_name='ga4gh.GetReferenceSetRequest.reference_set_id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=605,
  serialized_end=655,
)


_SEARCHREFERENCESREQUEST = _descriptor.Descriptor(
  name='SearchReferencesRequest',
  full_name='ga4gh.SearchReferencesRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='md5checksums', full_name='ga4gh.SearchReferencesRequest.md5checksums', index=0,
      number=1, type=9, cpp_type=9, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='accessions', full_name='ga4gh.SearchReferencesRequest.accessions', index=1,
      number=2, type=9, cpp_type=9, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='reference_set_id', full_name='ga4gh.SearchReferencesRequest.reference_set_id', index=2,
      number=3, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_token', full_name='ga4gh.SearchReferencesRequest.page_token', index=3,
      number=4, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_size', full_name='ga4gh.SearchReferencesRequest.page_size', index=4,
      number=5, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=658,
  serialized_end=790,
)


_SEARCHREFERENCESRESPONSE = _descriptor.Descriptor(
  name='SearchReferencesResponse',
  full_name='ga4gh.SearchReferencesResponse',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='references', full_name='ga4gh.SearchReferencesResponse.references', index=0,
      number=1, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='next_page_token', full_name='ga4gh.SearchReferencesResponse.next_page_token', index=1,
      number=2, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=792,
  serialized_end=881,
)


_GETREFERENCEREQUEST = _descriptor.Descriptor(
  name='GetReferenceRequest',
  full_name='ga4gh.GetReferenceRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='reference_id', full_name='ga4gh.GetReferenceRequest.reference_id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=883,
  serialized_end=926,
)


_LISTBASESREQUEST = _descriptor.Descriptor(
  name='ListBasesRequest',
  full_name='ga4gh.ListBasesRequest',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='reference_id', full_name='ga4gh.ListBasesRequest.reference_id', index=0,
      number=1, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='start', full_name='ga4gh.ListBasesRequest.start', index=1,
      number=2, type=3, cpp_type=2, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='end', full_name='ga4gh.ListBasesRequest.end', index=2,
      number=3, type=3, cpp_type=2, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_token', full_name='ga4gh.ListBasesRequest.page_token', index=3,
      number=4, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='page_size', full_name='ga4gh.ListBasesRequest.page_size', index=4,
      number=5, type=5, cpp_type=1, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=928,
  serialized_end=1035,
)


_LISTBASESRESPONSE = _descriptor.Descriptor(
  name='ListBasesResponse',
  full_name='ga4gh.ListBasesResponse',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  fields=[
    _descriptor.FieldDescriptor(
      name='offset', full_name='ga4gh.ListBasesResponse.offset', index=0,
      number=1, type=3, cpp_type=2, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='sequence', full_name='ga4gh.ListBasesResponse.sequence', index=1,
      number=2, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
    _descriptor.FieldDescriptor(
      name='next_page_token', full_name='ga4gh.ListBasesResponse.next_page_token', index=2,
      number=3, type=9, cpp_type=9, label=1,
      has_default_value=False, default_value=_b("").decode('utf-8'),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      options=None),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=1037,
  serialized_end=1115,
)

_SEARCHREFERENCESETSRESPONSE.fields_by_name['reference_sets'].message_type = _REFERENCESET
_SEARCHREFERENCESRESPONSE.fields_by_name['references'].message_type = _REFERENCE
DESCRIPTOR.message_types_by_name['Reference'] = _REFERENCE
DESCRIPTOR.message_types_by_name['ReferenceSet'] = _REFERENCESET
DESCRIPTOR.message_types_by_name['SearchReferenceSetsRequest'] = _SEARCHREFERENCESETSREQUEST
DESCRIPTOR.message_types_by_name['SearchReferenceSetsResponse'] = _SEARCHREFERENCESETSRESPONSE
DESCRIPTOR.message_types_by_name['GetReferenceSetRequest'] = _GETREFERENCESETREQUEST
DESCRIPTOR.message_types_by_name['SearchReferencesRequest'] = _SEARCHREFERENCESREQUEST
DESCRIPTOR.message_types_by_name['SearchReferencesResponse'] = _SEARCHREFERENCESRESPONSE
DESCRIPTOR.message_types_by_name['GetReferenceRequest'] = _GETREFERENCEREQUEST
DESCRIPTOR.message_types_by_name['ListBasesRequest'] = _LISTBASESREQUEST
DESCRIPTOR.message_types_by_name['ListBasesResponse'] = _LISTBASESRESPONSE

Reference = _reflection.GeneratedProtocolMessageType('Reference', (_message.Message,), dict(
  DESCRIPTOR = _REFERENCE,
  __module__ = 'proto.ga4gh.references_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.Reference)
  ))
_sym_db.RegisterMessage(Reference)

ReferenceSet = _reflection.GeneratedProtocolMessageType('ReferenceSet', (_message.Message,), dict(
  DESCRIPTOR = _REFERENCESET,
  __module__ = 'proto.ga4gh.references_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.ReferenceSet)
  ))
_sym_db.RegisterMessage(ReferenceSet)

SearchReferenceSetsRequest = _reflection.GeneratedProtocolMessageType('SearchReferenceSetsRequest', (_message.Message,), dict(
  DESCRIPTOR = _SEARCHREFERENCESETSREQUEST,
  __module__ = 'proto.ga4gh.references_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.SearchReferenceSetsRequest)
  ))
_sym_db.RegisterMessage(SearchReferenceSetsRequest)

SearchReferenceSetsResponse = _reflection.GeneratedProtocolMessageType('SearchReferenceSetsResponse', (_message.Message,), dict(
  DESCRIPTOR = _SEARCHREFERENCESETSRESPONSE,
  __module__ = 'proto.ga4gh.references_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.SearchReferenceSetsResponse)
  ))
_sym_db.RegisterMessage(SearchReferenceSetsResponse)

GetReferenceSetRequest = _reflection.GeneratedProtocolMessageType('GetReferenceSetRequest', (_message.Message,), dict(
  DESCRIPTOR = _GETREFERENCESETREQUEST,
  __module__ = 'proto.ga4gh.references_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.GetReferenceSetRequest)
  ))
_sym_db.RegisterMessage(GetReferenceSetRequest)

SearchReferencesRequest = _reflection.GeneratedProtocolMessageType('SearchReferencesRequest', (_message.Message,), dict(
  DESCRIPTOR = _SEARCHREFERENCESREQUEST,
  __module__ = 'proto.ga4gh.references_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.SearchReferencesRequest)
  ))
_sym_db.RegisterMessage(SearchReferencesRequest)

SearchReferencesResponse = _reflection.GeneratedProtocolMessageType('SearchReferencesResponse', (_message.Message,), dict(
  DESCRIPTOR = _SEARCHREFERENCESRESPONSE,
  __module__ = 'proto.ga4gh.references_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.SearchReferencesResponse)
  ))
_sym_db.RegisterMessage(SearchReferencesResponse)

GetReferenceRequest = _reflection.GeneratedProtocolMessageType('GetReferenceRequest', (_message.Message,), dict(
  DESCRIPTOR = _GETREFERENCEREQUEST,
  __module__ = 'proto.ga4gh.references_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.GetReferenceRequest)
  ))
_sym_db.RegisterMessage(GetReferenceRequest)

ListBasesRequest = _reflection.GeneratedProtocolMessageType('ListBasesRequest', (_message.Message,), dict(
  DESCRIPTOR = _LISTBASESREQUEST,
  __module__ = 'proto.ga4gh.references_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.ListBasesRequest)
  ))
_sym_db.RegisterMessage(ListBasesRequest)

ListBasesResponse = _reflection.GeneratedProtocolMessageType('ListBasesResponse', (_message.Message,), dict(
  DESCRIPTOR = _LISTBASESRESPONSE,
  __module__ = 'proto.ga4gh.references_pb2'
  # @@protoc_insertion_point(class_scope:ga4gh.ListBasesResponse)
  ))
_sym_db.RegisterMessage(ListBasesResponse)


# @@protoc_insertion_point(module_scope)
