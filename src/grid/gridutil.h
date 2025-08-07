//
/////
///// Class for managing MODFLOW-style arrays.  The class was written to
///// conserve memory so that CONSTANT arrays do not need to be allocated
///// with the full size.
/////
//
//template<class T>
//class ModflowArray
//{
//
//	/*!
//        \brief Create an instance of a ModflowArray object.
//
//         \param data The user provided input.  Can be a scalar value, which
//            is assigned to the entire array, a numpy array, or a list.
//
//         \param arrayshape The intended output shape of the array.
//         \param dtype The data type for the array (e.g. 'int', 'float', etc.).
//    */
//
//    ModflowArray()
//    {
//
//    }
//
//    /*
//    def __getitem__(self, key):
//        if isinstance(this->data, numpy.ndarray):
//            return this->data[key]
//        else:
//            return this->data
//    */
//
//	/// Convert the data to an array.
//    void as_ndarray(T * data, const list<int>& arrayshape)
//    {
//
//        //data is a numpy array
//        if( isinstance(this->data, numpy.ndarray) )
//        {
//
//            //if the shape is right, then just return the data
//            if( this->data.shape == this->arrayshape )
//                return this->data;
//
//            //if the array is 3D but only a 1D array was entered, then
//            //assume that there is one value for each layer
//            else if( len(this->data.shape) == 1 && arrayshape.size() == 3 )
//            {
//                nda = numpy.empty(this->arrayshape, dtype=this->dtype);
//                nlay = this->arrayshape[0];
//                for( k in range(nlay))
//                    nda[k, :, :] = this->data[k]
//
//                return nda;
//            }
//            //cannot convert array to an ndarray
//            else
//            {
//            	//we should create an exception here
//                cerr<<"Problem with converting data to a numpy array"<<endl;
//                assert(false);
//            }
//        }
//        //data is a list
//        else if( isinstance(this->data, list) )
//        {
//
//            //for a 1D array, just convert the list
//            if (len(this->arrayshape) == 1 &&
//               len(this->data) == this->arrayshape[0]):
//                nda = numpy.array(this->data)
//                return nda
//
//            else if( len(this->data) == this->arrayshape[0] )
//            {
//                nda = numpy.empty(this->arrayshape, dtype=this->dtype)
//                nlay = this->arrayshape[0]
//                for k in range(nlay):
//                    nda[k, :, :] = this->data[k]
//                return nda
//            }
//            else
//            {
//            	//we should create an exception here
//                cerr<<"Could not convert the input data list to a numpy array"<<endl;
//            	//we should create an exception here                assert(false);
//            }
//		}
//        //data is a scalar value
//        else
//        {
//        	//create a list of ones
//            nda = numpy.ones(this->arrayshape, dtype=this->dtype);
//
//            //copy the data from "data" to nda..
//            nda = nda * this->data;
//            return nda;
//        }
//
//}
//
///*
//class ModflowGridIndices()
//{
//    '''
//    Collection of methods that can be used to find cell indices for a
//    structured, but irregularly spaced MODFLOW grid.
//    '''
//
//    @staticmethod
//    def find_position_in_array(arr, x):
//        '''
//        If arr has x positions for the left edge of a cell, then return the cell
//        index containing x.
//
//        Arguments:
//
//            *arr*: A one dimensional array (such as Xe) that contains coordinates
//                for the left cell edge.
//
//            *x*: The x position to find in arr.
//        '''
//        jpos=None
//
//        if x == arr[-1]:
//            return len(arr) - 2
//
//        if x < min(arr[0], arr[-1]):
//            return None
//
//        if x > max(arr[0], arr[-1]):
//            return None
//
//        #go through each position
//        for j in range(len(arr)-1):
//            xl = arr[j]
//            xr = arr[j+1]
//            frac = (x - xl) / (xr - xl)
//            if 0. <= frac < 1.0:
//            #if min(xl, xr) <= x < max(xl, xr):
//                jpos = j
//                return jpos
//
//        return jpos
//
//    @staticmethod
//    def kij_from_nodenumber(nodenumber, nlay, nrow, ncol):
//        '''
//        Convert the modflow node number to a zero-based layer, row and column
//        format.  Return (k0, i0, j0).
//
//        Arguments:
//
//            *nodenumber*: The cell nodenumber, ranging from 1 to number of
//                nodes.
//
//            *nlay*: The number of layers.
//
//            *nrow*: The number of rows.
//
//            *ncol*: The number of columns.
//
//        '''
//        if nodenumber > nlay * nrow * ncol:
//            raise Exception('Error in function kij_from_nodenumber...')
//        n = nodenumber - 1
//        k = int( n / nrow / ncol )
//        i = int( (n - k * nrow * ncol) / ncol )
//        j = n - k * nrow * ncol - i * ncol
//        return (k, i, j)
//
//    @staticmethod
//    def nodenumber_from_kij(k, i, j, nrow, ncol):
//        '''
//        Calculate the nodenumber using the zero-based layer, row, and column
//        values.  The first node has a value of 1.
//
//        Arguments:
//
//            *k*: The model layer number as a zero-based value.
//
//            *i*: The model row number as a zero-based value.
//
//            *j*: The model column number as a zero-based value.
//
//            *nrow*: The number of model rows.
//
//            *ncol*: The number of model columns.
//        '''
//        return k * nrow * ncol + i * ncol + j + 1
//
//    @staticmethod
//    def nn0_from_kij(k, i, j, nrow, ncol):
//        '''
//        Calculate the zero-based nodenumber using the zero-based layer, row,
//         and column values.  The first node has a value of 0.
//
//        Arguments:
//
//            *k*: The model layer number as a zero-based value.
//
//            *i*: The model row number as a zero-based value.
//
//            *j*: The model column number as a zero-based value.
//
//            *nrow*: The number of model rows.
//
//            *ncol*: The number of model columns.
//        '''
//        return k * nrow * ncol + i * ncol + j
//
//    @staticmethod
//    def kij_from_nn0(n, nlay, nrow, ncol):
//        '''
//        Convert the node number to a zero-based layer, row and column
//        format.  Return (k0, i0, j0).
//
//        Arguments:
//
//            *nodenumber*: The cell nodenumber, ranging from 0 to number of
//                nodes - 1.
//
//            *nlay*: The number of layers.
//
//            *nrow*: The number of rows.
//
//            *ncol*: The number of columns.
//
//        '''
//        if n > nlay * nrow * ncol:
//            raise Exception('Error in function kij_from_nodenumber...')
//        k = int( n / nrow / ncol )
//        i = int( (n - k * nrow * ncol) / ncol )
//        j = n - k * nrow * ncol - i * ncol
//        return (k, i, j)
//}
//
//*/
