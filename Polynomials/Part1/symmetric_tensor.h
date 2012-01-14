#ifndef SYMMETRIC_TENSOR_H
#define SYMMETRIC_TENSOR_H

#include <array>

//Combinatorial helper methods
inline int factorial(int n) {
    int r = 1;
    for(int i=2; i<=n; ++i)
        r *= i;
    return r;
}

//Index iterator for tensors
template<int R, int D>
struct SymmetricIndexIterator {

    static_assert(R >= 0, "Rank must be nonnegative");
    static_assert(D > 0, "Dimension must be greater than 0");

    //Constants
    enum {
        rank = R,
        dimension = D,
    };

    //Type macros
    typedef std::array<int, dimension>  DegreeIndex;
    typedef std::array<int, rank>       TensorIndex;
    
    //Constructors
    SymmetricIndexIterator() : index_(0) {
        for(int i=0; i<rank; ++i)
            seq_[i] = 0;
    }
    SymmetricIndexIterator(int pos) : index_(pos) { }
    SymmetricIndexIterator(DegreeIndex degrees) :
        seq_(degree_to_tensor(degrees)),
        index_(degree_to_array(degrees)) {}
    
    //Sequence comparison
    inline bool operator==(SymmetricIndexIterator const& other) const  { return index_ == other.index_; }
    inline bool operator!=(SymmetricIndexIterator const& other) const  { return index_ != other.index_; }
    inline bool operator<(SymmetricIndexIterator const& other) const   { return index_ <  other.index_; }
    inline bool operator<=(SymmetricIndexIterator const& other) const  { return index_ <= other.index_; }
    inline bool operator>(SymmetricIndexIterator const& other) const   { return index_ >  other.index_; }
    inline bool operator>=(SymmetricIndexIterator const& other) const  { return index_ >= other.index_; }
    
    //Iteration
    SymmetricIndexIterator& operator++() {
        for(int i=0; i<rank; ++i) {
            if(seq_[i] < dimension - 1) {
                ++seq_[i];
                ++index_;
                for(int j=i-1; j>=0; --j) {
                    seq_[j] = seq_[i];
                }
                return *this;
            }
        }
        index_ = -1;
        return *this;
    }
    
    //Postfix iteration
    inline SymmetricIndexIterator operator++(int) {
        auto tmp = *this;
        ++*this;
        return tmp;
    }
    
    //Accessors
    inline DegreeIndex        degree_index() const { return tensor_to_degree(seq_); }
    inline TensorIndex const& tensor_index() const { return seq_; }    
    inline int                array_index()  const { return index_; }
    
    //Helper methods
    static inline int degree_to_array(DegreeIndex const& degrees) {
        int index = 0, sum = 0;
        for(int i=dimension-1; i>0; --i) {
            sum += degrees[i];
            int w = 1;
            for(int j=0; j<dimension-i; ++j) {
                w *= sum + j;
            }
            index += w / factorial(dimension - i);
        }
        return index;
    }

    static inline TensorIndex degree_to_tensor(DegreeIndex const& degrees) {
        TensorIndex result;
        int j = 0;
        for(int i=0; i<dimension; ++i) {
            for(int k=0; k<degrees[i]; ++k) {
                result[j++] = i;
            }
        }
        assert(j == rank);
        return result;
    }

    static inline DegreeIndex tensor_to_degree(TensorIndex const& seq) {
        DegreeIndex result;
        for(int i=0; i<dimension; ++i) {
            result[i] = 0;
        }
        for(int i=0; i<rank; ++i) {
            ++result[seq[i]];
        }
        return result;
    }

private:
    TensorIndex seq_;
    int index_;
};

#endif

