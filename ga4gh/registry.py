"""
The registry database for all the top-level objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import random
import hashlib
import datetime

import sqlalchemy
import sqlalchemy.exc
import sqlalchemy.orm as orm
import sqlalchemy.ext.declarative

import ga4gh.protocol as protocol
import ga4gh.exceptions as exceptions


def get_external_id():
    # TODO make this use SystemRandom and make a table specific prefix.
    return random.randint(0, 2**31)


def create_id_column():
    return sqlalchemy.Column(
        sqlalchemy.Integer, primary_key=True, autoincrement=False,
        default=get_external_id)


def create_timestamp_column():
    # TODO we should add a update default here as well.
    return sqlalchemy.Column(
        sqlalchemy.DateTime, nullable=False,
        default=datetime.datetime.now())


SqlAlchemyBase = sqlalchemy.ext.declarative.declarative_base()

class Info(SqlAlchemyBase):
    """
    A generic key-value table used to store info maps used in various
    places in the schema.
    """
    # TODO complete this. We probably need a more complicated
    # relationship because we can multiple values per key.
    __tablename__ = 'Info'

    id = sqlalchemy.Column(
        sqlalchemy.Integer, primary_key=True, autoincrement=True)
    key = sqlalchemy.Column(sqlalchemy.String, nullable=False, unique=True)
    value = sqlalchemy.Column(sqlalchemy.String, nullable=False)



class Dataset(SqlAlchemyBase):
    __tablename__ = 'Dataset'

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False, unique=True)
    description = sqlalchemy.Column(sqlalchemy.String, nullable=False)

    def __init__(self, name):
        self.name = name
        self.description = ""

    def get_protobuf(self):
        ret = protocol.Dataset()
        ret.id = str(self.id)
        ret.name = self.name
        ret.description = self.description
        return ret


class Reference(SqlAlchemyBase):
    __tablename__ = 'Reference'

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    md5checksum = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    length = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)
    source_uri = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    is_derived = sqlalchemy.Column(sqlalchemy.Boolean, nullable=False)
    source_divergence = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)

    reference_set_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("ReferenceSet.id"),
        nullable=False)
    reference_set = orm.relationship("ReferenceSet", back_populates="references")

    type = sqlalchemy.Column(sqlalchemy.String)
    __mapper_args__ = {
         'polymorphic_identity':'Reference',
         'polymorphic_on':type
    }
    __table_args__ = (
        # Reference names must be unique within a reference set.
        sqlalchemy.UniqueConstraint("reference_set_id", "name"),
    )

    def __init__(self, name):
        self.name = name
        self.source_divergence = 0
        self.source_uri = ""
        self.is_derived = False

    def get_protobuf(self):
        ret = protocol.Reference()
        ret.id = str(self.id)
        ret.name = self.name
        ret.md5checksum = self.md5checksum
        ret.length = self.length
        # We derive the ncbi_taxon_id from the reference set, since it
        # seems hard to imagine having a reference set containing references
        # from multiple species.
        ret.ncbi_taxon_id = self.reference_set.ncbi_taxon_id
        # ret.source_accessions.extend(self.source_accessions)
        ret.source_divergence = self.source_divergence
        ret.source_uri = self.source_uri
        return ret

class ReferenceSet(SqlAlchemyBase):
    __tablename__ = 'ReferenceSet'

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False, unique=True)
    description = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    assembly_id = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    md5checksum = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    is_derived = sqlalchemy.Column(sqlalchemy.Boolean, nullable=False)
    ncbi_taxon_id = sqlalchemy.Column(sqlalchemy.Integer, nullable=False)
    source_uri = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    # TODO source_acessions; should this be another table?
    references = orm.relationship(
        "Reference", back_populates="reference_set",
        cascade="all, delete, delete-orphan")
    type = sqlalchemy.Column(sqlalchemy.String)
    __mapper_args__ = {
         'polymorphic_identity':'ReferenceSet',
         'polymorphic_on':type
    }
    def __init__(self, name):
        self.name = name
        self.assembly_id = ""
        self.description = ""
        self.is_derived = False
        self.ncbi_taxon_id = 0
        self.source_uri = ""

    def _compute_md5checksum(self):
        """
        Computes the MD5 checksum for this reference set. This checksum is
        calculated by making a list of `Reference.md5checksum` for all
        `Reference`s in this set. We then sort this list, and take the
        MD5 hash of all the strings concatenated together.
        """
        checksums = sorted([ref.md5checksum for ref in self.references])
        self.md5checksum = hashlib.md5("".join(checksums)).hexdigest()

    def get_protobuf(self):
        ret = protocol.ReferenceSet()
        ret.id = str(self.id)
        ret.name = self.name
        ret.assembly_id = self.assembly_id
        ret.description = self.description
        ret.is_derived = self.is_derived
        ret.md5checksum = self.md5checksum
        ret.ncbi_taxon_id = self.ncbi_taxon_id
        # ret.source_accessions.extend(self.getSourceAccessions())
        ret.source_uri = self.source_uri
        return ret

class VariantSetMetadata(SqlAlchemyBase):
    __tablename__ = "VariantSetMetadata"

    id = create_id_column()
    key = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    value = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    type = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    number = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    description = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    variant_set_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("VariantSet.id"),
        nullable=False)
    variant_set = orm.relationship(
        "VariantSet", back_populates="variant_set_metadata")

    def __init__(self, key=None, value="", type="", number="", description=""):
        self.key = key
        self.value = value
        self.type = type
        self.number = number
        self.description = description

    def get_protobuf(self):
        ret = protocol.VariantSetMetadata()
        ret.id = str(self.id)
        ret.key = self.key
        ret.value = self.value
        ret.type = self.type
        ret.number = self.number
        ret.description = self.description
        return ret

# There is a many-to-many association between VariantSets and CallSets,
# so we need an association table.
call_set_association_table = sqlalchemy.Table(
    'call_set_association', SqlAlchemyBase.metadata,
    sqlalchemy.Column(
        'variant_set_id', sqlalchemy.Integer,
        sqlalchemy.ForeignKey('VariantSet.id')),
    sqlalchemy.Column(
        'call_set_id', sqlalchemy.Integer,
        sqlalchemy.ForeignKey('CallSet.id')))


class CallSet(SqlAlchemyBase):
    __tablename__ = 'CallSet'

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    sample_id = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    creation_timestamp = create_timestamp_column()
    update_timestamp = create_timestamp_column()
    variant_sets = orm.relationship(
        "VariantSet", secondary=call_set_association_table,
        back_populates="call_sets")

    def get_protobuf(self):
        ret = protocol.CallSet()
        ret.id = str(self.id)
        ret.name = self.name
        ret.sample_id = self.sample_id
        ret.created = protocol.datetime_to_milliseconds(
            self.creation_timestamp)
        ret.updated = protocol.datetime_to_milliseconds(
            self.update_timestamp)
        ret.variant_set_ids.extend(str(vs.id) for vs in self.variant_sets)
        # TODO create a generic Info table that can hold these values
        # for a variety of different tables. Or it may be simpler to create
        # a per-class table.
        # for key in self._info:
        #     gaCallSet.info[key].values.extend(_encodeValue(self._info[key]))
        return ret

    def __init__(self, name=None, sample_id=None):
        self.name = name
        self.sample_id = sample_id


class VariantSet(SqlAlchemyBase):
    __tablename__ = 'VariantSet'

    id = create_id_column()
    name = sqlalchemy.Column(sqlalchemy.String, nullable=False)
    creation_timestamp = create_timestamp_column()
    update_timestamp = create_timestamp_column()

    dataset_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("Dataset.id"),
        nullable=False)
    # dataset = orm.relationship(
    #     "Dataset", back_populates="dataset", single_parent=True,
    #     cascade="all, delete, delete-orphan")
    reference_set_id = sqlalchemy.Column(
        sqlalchemy.Integer, sqlalchemy.ForeignKey("ReferenceSet.id"),
        nullable=False)
    # reference_set = orm.relationship(
    #     "ReferenceSet", back_populates="reference_set")

    variant_set_metadata = orm.relationship(
        "VariantSetMetadata", cascade="all, delete, delete-orphan")
    call_sets = orm.relationship(
        "CallSet", secondary=call_set_association_table,
        back_populates="variant_sets", lazy="dynamic")

    type = sqlalchemy.Column(sqlalchemy.String)
    __mapper_args__ = {
         'polymorphic_identity':'VariantSet',
         'polymorphic_on':type
    }
    __table_args__ = (
        # VariantSet names must be unique within a dataset
        sqlalchemy.UniqueConstraint("dataset_id", "name"),
    )

    def __init__(self, name):
        self.name = name

    def get_protobuf(self):
        ret = protocol.VariantSet()
        ret.id = str(self.id)
        ret.name = self.name
        # ret.created = self.creation_timestamp
        # ret.updated = self.update_timestamp
        ret.dataset_id = str(self.dataset_id)
        ret.reference_set_id = str(self.reference_set_id)
        metadata = [m.get_protobuf() for m in self.variant_set_metadata]
        ret.metadata.extend(metadata)
        return ret



class RegistryDb(object):
    """
    The database representing all of the top-level objects in a given
    data repository.
    """
    def __init__(self, db_url):
        self._db_url = db_url
        self._engine = None
        self._session = None

    def get_new_id(self):
        # TODO use an instance of SystemRandom to do this and check against the
        # global IDs table. We should not have the possibility of having two
        # objects with the same ID. We could avoid the problem by using the
        # first 4 bits to identify the type.
        return random.randint(0, 2**31)

    # Miscellaneous house-keeping methods for the Registry DB.

    def exists(self):
        """
        Returns True if there is a fully initialised registry in the DB url.
        """
        return True

    def open(self):
        """
        Opens this registry DB.
        """
        self._engine = sqlalchemy.create_engine(self._db_url, echo=False)
        Session = orm.sessionmaker()
        Session.configure(bind=self._engine)
        self._session = Session()

    def commit(self):
        """
        Commits any changes made to the repo. It is an error to call
        this function if the repo is not opened in write-mode.
        """
        self._session.commit()

    def close(self):
        """
        Closes this repo.
        """
        self._session.close()
        self._session = None
        self._engine = None

    def delete(self):
        """
        Drops all tables in the DB.
        """
        self.open()
        metadata = datamodel.SqlAlchemyBase.metadata
        metadata.drop_all(self._engine)
        self.close()

    def initialise(self):
        """
        Initialise this data repostitory by creating the registry DB.
        """
        metadata = SqlAlchemyBase.metadata
        metadata.drop_all(self._engine)
        metadata.create_all(self._engine)


    # Top level API to add objects to the data repository.

    def add_reference_set(self, reference_set):
        """
        Adds the specified reference set to this data repository.
        """
        try:
            self._session.add(reference_set)
            self._session.commit()
        except sqlalchemy.exc.IntegrityError as ie:
            raise exceptions.DuplicateNameException(reference_set.name)

    def add_dataset(self, dataset):
        """
        Adds the specified dataset to this data repository.
        """
        try:
            self._session.add(dataset)
            self._session.commit()
        except sqlalchemy.exc.IntegrityError as ie:
            raise exceptions.DuplicateNameException(dataset.name)

    def add_variant_set(self, variant_set):
        """
        Adds the specified variant_set to this data repository.
        """
        try:
            self._session.add(variant_set)
            self._session.commit()
        except sqlalchemy.exc.IntegrityError as ie:
            raise exceptions.DuplicateNameException(variant_set.name)

    # Object accessors by name.

    def get_dataset_by_name(self, name):
        """
        Returns the dataset with the specified name.
        """
        result = self._session.query(Dataset).filter(
            Dataset.name == name).first()
        if result is None:
            raise exceptions.DatasetNameNotFoundException(name)
        return result

    def get_reference_set_by_name(self, name):
        """
        Returns the reference_set with the specified name.
        """
        result = self._session.query(ReferenceSet).filter(
            ReferenceSet.name == name).first()
        if result is None:
            raise exceptions.ReferenceSetNameNotFoundException(name)
        return result

    # Object accessors by ID. These are run when the corresponding GET request
    # is received.

    def get_call_set(self, id_):
        """
        Retuns the CallSet with the specified ID, or raises a
        VariantSetNotFoundException if it does not exist.
        """
        result = self._session.query(CallSet).filter(CallSet.id == id_).first()
        if result is None:
            raise exceptions.CallSetNotFoundException(id_)
        return result

    def get_variant_set(self, id_):
        """
        Retuns the VariantSet with the specified ID, or raises a
        VariantSetNotFoundException if it does not exist.
        """
        result = self._session.query(VariantSet).filter(
            VariantSet.id == id_).first()
        if result is None:
            raise exceptions.VariantSetNotFoundException(id_)
        return result

    def get_reference_set(self, id_):
        """
        Retuns the ReferenceSet with the specified ID, or raises a
        ReferenceSetNotFoundException if it does not exist.
        """
        result = self._session.query(ReferenceSet).filter(
            ReferenceSet.id == id_).first()
        if result is None:
            raise exceptions.ReferenceSetNotFoundException(id_)
        return result

    def get_reference(self, id_):
        """
        Returns the Reference with the specified ID or raises a
        ReferenceNotFoundException if it does not exist.
        """
        result = self._session.query(Reference).filter(
            Reference.id == id_).first()
        if result is None:
            raise exceptions.ReferenceNotFoundException(id_)
        return result

    # General getters to iterate over all objects in a class.

    def get_reference_sets(self):
        return self._session.query(ReferenceSet).all()

    # Getters to provide queries corresponding to the external search requests.

    def get_references_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchReferencesRequest.
        """
        query = self._session.query(Reference).filter(
            Reference.reference_set_id == request.reference_set_id)
        if request.md5checksum:
            query = query.filter(Reference.md5checksum == request.md5checksum)
        return query

    def get_reference_sets_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchReferenceSetsRequest.
        """
        query = self._session.query(ReferenceSet)
        if request.md5checksum:
            query = query.filter(ReferenceSet.md5checksum == request.md5checksum)
        if request.assembly_id:
            query = query.filter(ReferenceSet.assembly_id == request.assembly_id)
        return query

    def get_datasets_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchDatasetsRequest.
        """
        query = self._session.query(Dataset)
        return query

    def get_variant_sets_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchVariantSets request.
        """
        query = self._session.query(VariantSet)
        return query

    def get_call_sets_search_query(self, request):
        """
        Returns the query object representing all the rows in the specified
        SearchCallSets request.
        """
        query = self._session.query(CallSet)
        query = query.filter(CallSet.variant_sets.any(id=request.variant_set_id))
        if request.name:
            query = query.filter(CallSet.name == request.name)
        return query
